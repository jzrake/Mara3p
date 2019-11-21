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
#include "model_atmospheres.hpp"
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
    .item("dfi",                  1.5)   // output frequency (constant multiplier)
    .item("rk_order",               2)   // Runge-Kutta order (1 or 2)
    .item("cfl",                  0.5)   // courant number
    .item("plm_theta",            1.5)   // PLM parameter
    .item("router",               1e3)   // outer boundary radius
    .item("move",                   1)   // whether to move the cells
    .item("ambient_medium_index", 2.0)   // index in A * (r / r0)^(-index)
    .item("ambient_medium_norm",  1.0)   //     A in A * (r / r0)^(-index)
    .item("Lwind",                1.0)   // wind luminosity at t=1
    .item("rwind",                1.0)   // radius out to which an initial wind profile is set
    .item("uwind",               10.0)   // radial gamma-beta of the wind
    .item("twind",                1.0)   // time scale the wind lasts for
    .item("envelop",                1)   // include a relativistic envelop as part of the ambient medium
    .item("delay",               10.0);  // the delay time [inner boundary light-crossing time] of the cloud
}

static const auto gamma_law_index   = 4. / 3;
static const auto temperature_floor = 0.0;




//=============================================================================
using solution_t            = mara::state_with_vertices_t<srhd::conserved_t>;
using solution_with_tasks_t = std::pair<solution_t, control::task_t>;
using timed_state_pair_t    = control::timed_pair_t<solution_with_tasks_t>;
using mara::spherical_mesh_geometry_t;

inline auto solution(solution_with_tasks_t p)
{
    return p.first;
}

inline auto tasks(solution_with_tasks_t p)
{
    return p.second;
}




//=============================================================================
auto cold_power_law_medium(const mara::config_t& run_config)
{
    return mara::cold_power_law_medium_t()
    .with_density_index  (run_config.get_double("ambient_medium_index"))
    .with_density_at_base(run_config.get_double("ambient_medium_norm"));
}

auto cold_wind(const mara::config_t& run_config)
{
    return mara::cold_wind_model_t()
    .with_kinetic_luminosity(run_config.get_double("Lwind"))
    .with_time_scale        (run_config.get_double("twind"))
    .with_gamma_beta        (run_config.get_double("uwind"))
    .with_solid_angle       (4 * M_PI);
}

auto cloud_with_envelop(const mara::config_t& run_config)
{
    return mara::cloud_with_envelop_model_t();
}




auto fast_rise_exponential_decay(dimensional::unit_time onset_time, dimensional::unit_time decay_time)
{
    return [=] (dimensional::unit_time t)
    {
        auto x = 2.0 * (t - onset_time) / decay_time;
        return 1.0 / (4 * std::exp(-2)) * x * x * std::exp(-x) * (x < 0.0 ? 0.0 : 1.0);
    };
}




auto envelop_mdot_index       = 4.0;
auto envelop_gamma_beta_index = 0.22; // psi parameter
auto envelop_nominal_time     = dimensional::unit_time(1.0);
auto envelop_nominal_mdot     = dimensional::unit_mass_rate(1.0);
auto envelop_nominal_mass     = dimensional::unit_mass(1.0);
auto envelop_gamma_beta       = dimensional::unit_scalar(20.0);

auto engine_gamma_beta     = dimensional::unit_scalar(20.0);
auto engine_onset_time     = dimensional::unit_time(50.0);
auto engine_duration       = dimensional::unit_time(100.0);




auto wind_mass_loss_rate(const mara::config_t& run_config)
{
    return [] (auto t)
    {
        // auto m0 = envelop_nominal_mass;
        auto t0 = envelop_nominal_time;
        auto al = envelop_mdot_index;
        return t < engine_onset_time
        ? dimensional::unit_mass_rate(1.0) * std::pow(t / t0, al)
        : dimensional::unit_mass_rate(1e4);
    };
}

auto wind_gamma_beta(const mara::config_t& run_config)
{
    auto integrated_envelop_mdot = [] (dimensional::unit_time t)
    {
        auto m0 = envelop_nominal_mass;
        auto md = envelop_nominal_mdot;
        auto t0 = envelop_nominal_time;
        auto al = envelop_mdot_index;
        return m0 + md * t * std::pow(t / t0, al) / (1 + al);
        // return m0 * (1.0 + std::erf(t / t0 - 1.0));
    };

    return [integrated_envelop_mdot] (auto t)
    {
        auto psi = envelop_gamma_beta_index;
        auto m  = integrated_envelop_mdot(t);
        auto m0 = integrated_envelop_mdot(envelop_nominal_time);
        auto u0 = envelop_gamma_beta;

        return t < engine_onset_time
        ? u0 * std::pow(m / m0, -psi)
        : engine_gamma_beta * fast_rise_exponential_decay(engine_onset_time, engine_duration)(t);
    };
}




//=============================================================================
auto wind_profile(const mara::config_t& run_config, dimensional::unit_length r, dimensional::unit_time t)
{
    return mara::time_varying_cold_wind_t()
    .with_mass_loss_rate(wind_mass_loss_rate(run_config))
    .with_gamma_beta    (wind_gamma_beta    (run_config))
    .primitive_srhd(r, t);
}

auto initial_condition(const mara::config_t& run_config)
{
    return [run_config] (dimensional::unit_length r) { return wind_profile(run_config, r, 1.0); };
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
        : std::make_pair(srhd::riemann_hllc(pl, pr, nhat, gamma_law_index), dimensional::unit_velocity(0.0));
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
    | nd::construct<dimensional::unit_length>();
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
solution_t remesh(const mara::config_t& run_config, solution_t solution)
{
    if (front(solution.vertices) > dimensional::unit_length(1.0 + 1.0 / run_config.get_int("nr")))
    {
        auto xv = nd::concat(nd::from(dimensional::unit_length(1.0)), solution.vertices, 0) | nd::to_shared();
        auto xc = spherical_mesh_geometry_t::cell_centers(xv);
        auto dv = spherical_mesh_geometry_t::cell_volumes(xv);
        auto bp = wind_profile(run_config, front(xc), solution.time);
        auto bu = srhd::conserved_density(bp, gamma_law_index) * front(dv);
        auto u1 = nd::concat(nd::from(bu), solution.conserved, 0) | nd::to_shared();

        return {
            solution.iteration,
            solution.time,
            xv,
            u1,
        };
    }
    return solution;
}

solution_t advance(const mara::config_t& run_config, solution_t solution)
{
    auto base = [&run_config, dt = time_step(run_config, solution)] (solution_t soln)
    {
        auto move_cells          = run_config.get_int("move");
        auto plm_theta           = run_config.get_double("plm_theta");
        auto xhat                = geometric::unit_vector_on_axis(1);
        auto inner_boundary_prim = wind_profile     (run_config, front(soln.vertices), soln.time);
        auto outer_boundary_prim = initial_condition(run_config)(back (soln.vertices));
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

control::task_t advance(const mara::config_t& run_config, control::task_t task, dimensional::unit_time time)
{
    auto reference_time = time < engine_onset_time ? dimensional::unit_time(0.0) : engine_onset_time;
    return jump(task, time, run_config.get_double("dfi"), reference_time);
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
        return solution(control::last_state(p)).time <= dimensional::unit_time(tfinal);
    };
}

void print_run_loop(const mara::config_t& run_config, timed_state_pair_t p)
{
    auto soln = solution(control::this_state(p));
    auto nz = size(soln.vertices);
    auto us = control::microseconds_separating(p);

    std::printf("[%07lu] t=%.3lf dt=%.2e zones=%lu Mzps=%.2lf ",
        long(soln.iteration),
        soln.time.value,
        time_step(run_config, soln).value,
        nz,
        nz / us);

    auto dr = soln.vertices | nd::adjacent_diff();

    std::printf("| min(dr)=%.2e @ %lu max(dr)=%.2e @ %lu\n",
        nd::min(dr).value, nd::argmin(dr)[0],
        nd::max(dr).value, nd::argmax(dr)[0]);
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
