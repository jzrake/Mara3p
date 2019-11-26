/**
 ==============================================================================
 Copyright 2019, Jonathan Zrake

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 ==============================================================================
*/




// #define SRHD_NO_EXCEPTIONS
#include "app_config.hpp"
#include "app_control.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_config.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_state_templates.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_sequence.hpp"
#include "core_util.hpp"
#include "model_wind.hpp"
#include "physics_srhd.hpp"
#include "scheme_mesh_geometry.hpp"
#include "scheme_moving_mesh.hpp"
#include "scheme_plm_gradient.hpp"




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("nr",                   256)   // number of radial zones, per decade
    .item("tfinal",              10.0)   // time to stop the simulation
    .item("print",                 10)   // the number of iterations between terminal outputs
    .item("dfi",                  1.5)   // output interval (constant multiplier)
    .item("rk_order",               2)   // Runge-Kutta order (1, 2, or 3)
    .item("cfl",                  0.5)   // courant number
    .item("mindr",               1e-3)   // minimum cell length to impose in remeshing
    .item("plm_theta",            1.5)   // PLM parameter
    .item("router",               1e3)   // outer boundary radius
    .item("move",                   1)   // whether to move the cells
    .item("envelop_mdot_index",   4.0)   // alpha in envelop mdot(t) = (t / t0)^alpha
    .item("envelop_u_index",     0.22)   // psi in envelop u(m) = u1 (m / m1)^(-psi)
    .item("envelop_u",          10.00)   // maximum gamma-beta in outer envelop
    .item("engine_mdot",          1e4)   // engine mass rate
    .item("engine_onset",        50.0)   // the engine onset time [inner boundary light-crossing time]
    .item("engine_duration",    100.0)   // the engine duration   [inner boundary light-crossing time]
    .item("engine_u",            10.0);  // the engine gamma-beta
}

static const auto gamma_law_index   = 4. / 3;
static const auto temperature_floor = 0.0;




//=============================================================================
using solution_t            = mara::state_with_vertices_t<srhd::conserved_t>;
using solution_with_tasks_t = std::pair<solution_t, control::task_t>;
using timed_state_pair_t    = control::timed_pair_t<solution_with_tasks_t>;
using mara::spherical_mesh_geometry_t;
using namespace dimensional;

inline auto solution(solution_with_tasks_t p)
{
    return p.first;
}

inline auto tasks(solution_with_tasks_t p)
{
    return p.second;
}




auto wind_mass_loss_rate(const mara::config_t& run_config)
{
    auto envelop_mdot = unit_mass_rate(1.0);
    auto engine_mdot  = unit_mass_rate(run_config.get_double("engine_mdot"));
    auto engine_onset = unit_time(run_config.get_double("engine_onset"));
    auto alpha        = run_config.get_double("envelop_mdot_index");

    return [envelop_mdot, engine_mdot, alpha, engine_onset] (auto t)
    {
        return t < engine_onset
        ? envelop_mdot * std::pow(t / unit_time(1.0), alpha)
        : engine_mdot;
    };
}

auto wind_gamma_beta(const mara::config_t& run_config)
{
    auto m0              = unit_mass(1.0);
    auto md              = unit_mass_rate(1.0);
    auto engine_onset    = unit_time(run_config.get_double("engine_onset"));
    auto engine_duration = unit_time(run_config.get_double("engine_duration"));
    auto engine_u        = unit_scalar(run_config.get_double("engine_u"));
    auto envelop_u       = unit_scalar(run_config.get_double("envelop_u"));
    auto fred            = mara::fast_rise_exponential_decay(engine_onset, engine_duration);
    auto psi             = run_config.get_double("envelop_u_index");
    auto alpha           = run_config.get_double("envelop_mdot_index");

    auto integrated_envelop_mdot = [m0, md, alpha] (unit_time t)
    {
        return m0 + md * t * std::pow(t / unit_time(1.0), alpha) / (1 + alpha);
    };

    return [integrated_envelop_mdot, envelop_u, engine_u, fred, psi] (auto t)
    {
        auto m0   = integrated_envelop_mdot(1.0);
        auto m    = integrated_envelop_mdot(t);
        return envelop_u * std::pow(m / m0, -psi) + engine_u * fred(t);
    };
}




//=============================================================================
auto wind_profile(const mara::config_t& run_config, unit_length r, unit_time t)
{
    return mara::cold_relativistic_wind_t()
    .with_mass_loss_rate(wind_mass_loss_rate(run_config))
    .with_gamma_beta    (wind_gamma_beta    (run_config))
    .primitive(r, t);
}

auto initial_condition(const mara::config_t& run_config)
{
    return [run_config] (unit_length r) { return wind_profile(run_config, r, 1.0); };
}

auto recover_primitive()
{
    return [] (auto u) { return srhd::recover_primitive(u, gamma_law_index, temperature_floor); };
}

auto riemann_solver_for(geometric::unit_vector_t nhat, bool move)
{
    auto contact_mode = srhd::riemann_solver_mode_hllc_fluxes_across_contact_t();

    return util::apply_to([nhat, move, contact_mode] (auto pl, auto pr)
    {
        return move
        ? srhd::riemann_solver(pl, pr, nhat, gamma_law_index, contact_mode)
        : std::make_pair(srhd::riemann_hllc(pl, pr, nhat, gamma_law_index), unit_velocity(0.0));
    });
}

auto time_step(const mara::config_t& run_config, solution_t state)
{
    auto cfl = run_config.get_double("cfl");
    return cfl * nd::min(state.vertices | nd::adjacent_diff()) / srhd::light_speed;
}




//=============================================================================
auto initial_vertices(const mara::config_t& run_config)
{
    auto router     = run_config.get_double("router");
    auto cell_count = run_config.get_int("nr") * int(std::log10(router));

    return nd::linspace(0.0, std::log10(router), cell_count + 1)
    | nd::map([] (auto log10r) { return std::pow(10.0, log10r); })
    | nd::construct<unit_length>();
}

solution_t initial_solution_state(const mara::config_t& run_config)
{
    auto xv = initial_vertices(run_config) | nd::to_shared();
    auto dv = spherical_mesh_geometry_t::cell_volumes(xv);
    auto p0 = xv | nd::adjacent_mean() | nd::map(initial_condition(run_config));
    auto u0 = p0 | nd::map([] (auto p) { return srhd::conserved_density(p, gamma_law_index); });
    return {0, 1.0, xv, (u0 * dv) | nd::to_shared()};
}

solution_with_tasks_t initial_app_state(const mara::config_t& run_config)
{
    return std::pair(initial_solution_state(run_config), control::task("write_diagnostics", 1.0));
}




//=============================================================================
solution_t add_inner_cell(const mara::config_t& run_config, solution_t solution)
{
    if (front(solution.vertices) > unit_length(1.0 + 1.0 / run_config.get_int("nr")))
    {
        auto x1 = nd::concat(nd::from(unit_length(1.0)), solution.vertices, 0) | nd::to_shared();
        auto xc = spherical_mesh_geometry_t::cell_centers(x1);
        auto dv = spherical_mesh_geometry_t::cell_volumes(x1);
        auto bp = wind_profile(run_config, front(xc), solution.time);
        auto bu = srhd::conserved_density(bp, gamma_law_index) * front(dv);
        auto u1 = nd::concat(nd::from(bu), solution.conserved, 0) | nd::to_shared();

        return {
            solution.iteration,
            solution.time,
            x1,
            u1,
        };
    }
    return solution;
}

solution_t remove_smallest_cell(const mara::config_t& run_config, solution_t solution)
{
    auto mindr = unit_length(run_config.get_double("mindr"));
    auto dr    = spherical_mesh_geometry_t::cell_spacings(solution.vertices);
    auto imin  = nd::argmin(dr)[0];

    if (imin > 0 && imin < size(solution.conserved) && dr(imin) < mindr)
    {
        auto [x1, u1] = nd::remove_partition(solution.vertices, solution.conserved, imin);

        return {
            solution.iteration,
            solution.time,
            x1 | nd::to_shared(),
            u1 | nd::to_shared(),
        };
    }
    return solution;
}

solution_t remesh(const mara::config_t& run_config, solution_t solution)
{
    solution = add_inner_cell      (run_config, solution);
    solution = remove_smallest_cell(run_config, solution);
    return solution;
}

solution_t advance(const mara::config_t& run_config, solution_t solution)
{
    auto base = [&run_config, dt = time_step(run_config, solution)] (solution_t soln)
    {
        auto move_cells          = run_config.get_int("move");
        auto plm_theta           = run_config.get_double("plm_theta");
        auto xhat                = geometric::unit_vector_on_axis(1);
        auto inner_boundary_prim = wind_profile(run_config, front(soln.vertices), soln.time);
        auto outer_boundary_prim = wind_profile(run_config, back (soln.vertices), 1.0);
        auto mesh_geometry       = spherical_mesh_geometry_t();
        auto source_terms        = util::apply_to([] (auto p, auto x)
        {
            return srhd::spherical_geometry_source_terms(p, x, M_PI / 2, gamma_law_index);
        });
        return mara::advance(
            soln,
            dt,
            inner_boundary_prim,
            outer_boundary_prim,
            riemann_solver_for(xhat, move_cells),
            recover_primitive(),
            source_terms,
            mesh_geometry,
            plm_theta);
    };
    return remesh(run_config, control::advance_runge_kutta(base, run_config.get_int("rk_order"), solution));
}

control::task_t advance(const mara::config_t& run_config, control::task_t task, unit_time time)
{
    auto engine_onset = unit_time(run_config.get_double("engine_onset"));
    auto ref_time = time < engine_onset ? unit_time(0.0) : engine_onset;
    return jump(task, time, run_config.get_double("dfi"), ref_time);
}

auto advance(const mara::config_t& run_config)
{
    return util::apply_to([run_config] (solution_t state, control::task_t task)
    {
        return std::pair(
            advance(run_config, state),
            advance(run_config, task, state.time));
    });
}




//=============================================================================
void write_diagnostics(const mara::config_t& run_config, solution_t state, unsigned long count)
{
    auto fname     = util::format("sedov.%04lu.h5", count);
    auto dv        = spherical_mesh_geometry_t::cell_volumes(state.vertices);
    auto prim      = state.conserved | nd::divide(dv) | nd::map(recover_primitive()) | nd::to_shared();
    auto file      = h5::File(fname, "w");

    std::printf("Write diagnostics %s\n", fname.data());

    h5::write(file, "vertices", state.vertices);
    h5::write(file, "primitive", prim);
    h5::write(file, "time", state.time);
    h5::write(file, "run_config", run_config);
}

auto should_continue(const mara::config_t& run_config)
{
    return [tfinal = run_config.get_double("tfinal")] (timed_state_pair_t p)
    {
        return solution(control::last_state(p)).time <= unit_time(tfinal);
    };
}

void print_run_loop(const mara::config_t& run_config, timed_state_pair_t p)
{
    auto soln = solution(control::this_state(p));
    auto nz = size(soln.vertices);
    auto us = control::microseconds_separating(p);

    std::printf("[%07lu] t=%.3lf dt=%.2e zones=%lu Mzps=%.2lf\n",
        long(soln.iteration),
        soln.time.value,
        time_step(run_config, soln).value,
        nz,
        nz / us);

    // auto dr = soln.vertices | nd::adjacent_diff();
    // std::printf("| min(dr)=%.2e @ %lu max(dr)=%.2e @ %lu\n",
    //     nd::min(dr).value, nd::argmin(dr)[0],
    //     nd::max(dr).value, nd::argmax(dr)[0]);
}

auto side_effects(const mara::config_t& run_config, timed_state_pair_t p)
{
    auto this_soln = solution(control::this_state(p));
    auto last_soln = solution(control::last_state(p));
    auto this_task = tasks(control::this_state(p));
    auto last_task = tasks(control::last_state(p));

    if (this_task.count != last_task.count)
        write_diagnostics(run_config, last_soln, last_task.count);

    if (long(this_soln.iteration) % run_config.get_int("print") == 0)
        print_run_loop(run_config, p);            
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto run_config = config_template()
    .create()
    .update(mara::argv_to_string_map(argc, argv));

    mara::pretty_print(std::cout, "config", run_config);

    auto simulation = seq::generate(initial_app_state(run_config), advance(run_config))
    | seq::pair_with(control::time_point_sequence())
    | seq::window()
    | seq::take_while(should_continue(run_config));

    for (auto state_pair : simulation)
    {
        side_effects(run_config, state_pair);
    }
    return 0;
}
