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




#include "app_config.hpp"
#include "app_control.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_state_templates.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_sequence.hpp"
#include "physics_srhd.hpp"
#include "scheme_plm_gradient.hpp"




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("nr",                 256)   // number of radial zones, per decade
    .item("tfinal",             1.0)   // time to stop the simulation
    .item("cfl",                0.5)   // courant number
    .item("outer_radius",     100.0);  // outer boundary radius
}

static const auto gamma_law_index = 4. / 3;




//=============================================================================
using state_t = mara::state_with_vertices_t<srhd::conserved_t>;
using timed_state_pair_t = control::timed_pair_t<state_t>;




//=============================================================================
template<typename FunctionType>
auto apply_to(FunctionType function)
{
    return [function] (auto t) { return std::apply(function, t); };
}




//=============================================================================
auto initial_vertices(nd::uint cell_count)
{
    return nd::linspace(0.0, 2.0, cell_count + 1)
    | nd::map([] (auto log10r) { return std::pow(10.0, log10r); })
    | nd::construct<dimensional::unit_length>();
}

auto shell_volume(dimensional::unit_length r0, dimensional::unit_length r1)
{
    return (r1 * r1 * r1 - r0 * r0 * r0) / 3.0;
}




//=============================================================================
auto cell_centers(state_t state)
{
    return state.vertices | nd::adjacent_mean();
}

auto cell_spacing(state_t state)
{
    return state.vertices | nd::adjacent_diff();
}

auto cell_volumes(state_t state)
{
    return state.vertices | nd::adjacent_zip() | nd::map(apply_to(shell_volume));
}

auto face_areas(state_t state)
{
    return state.vertices | nd::map([] (auto r) { return r * r; });
}

auto time_step(state_t state)
{
    return 0.33 * nd::min(state.vertices | nd::adjacent_diff()) / srhd::light_speed;
}




//=============================================================================
auto initial_condition(const mara::config_t& run_config, dimensional::unit_length r)
{
    auto r0 = dimensional::unit_length(1.0);
    auto d_wind = dimensional::unit_mass_density(1.0) * std::pow(r / r0, -2.0);
    auto u_wind = 2.0;
    auto p_wind = d_wind * 1e-3 * dimensional::pow<2>(srhd::light_speed);
    auto d_upst = 1e-3;
    auto p_upst = 1e-3 * d_upst;

    auto wind = srhd::primitive(d_wind, u_wind, 0.0, 0.0, p_wind);
    auto upst = srhd::primitive(d_upst,    0.0, 0.0, 0.0, p_upst);

    return r < dimensional::unit_length(10.0) ? wind : upst;
}

auto initial_condition(const mara::config_t& run_config)
{
    return [run_config] (dimensional::unit_length r)
    {
        return initial_condition(run_config, r);
    };
}

auto riemann_solver_for(geometric::unit_vector_t nhat)
{
    return [nhat] (auto t)
    {
        auto [pf, vf] = t;
        auto [pl, pr] = pf;
        return srhd::riemann_hlle(pl, pr, nhat, vf, gamma_law_index);
    };
}

auto recover_primitive()
{
    return [] (auto u) { return srhd::recover_primitive(u, gamma_law_index); };
}




//=============================================================================
state_t initial_state(const mara::config_t& run_config)
{
    auto xv = initial_vertices(run_config.get_int("nr") * int(std::log10(run_config.get_double("outer_radius"))));
    auto dv = xv | nd::adjacent_zip() | nd::map(apply_to(shell_volume));
    auto p0 = xv | nd::adjacent_mean() | nd::map(initial_condition(run_config));
    auto u0 = p0 | nd::map([] (auto p) { return srhd::conserved_density(p, gamma_law_index); });
    return {0.0, 0, xv | nd::to_shared(), (u0 * dv) | nd::to_shared()};
}




//=============================================================================
state_t advance(const mara::config_t& run_config, state_t state)
{
    /*
        p0  :=    x   x   x   x   x
        s0  :=    x   x   x   x   x
        gx  :=    o   x   x   x   o      (o = extend zeros)
        pf  :=      |   |   |   |
        vf  :=      |   |   |   |
        f0  :=  :   |   |   |   |   ;    (: = extend-uniform-lower ; = extend-zero-gradient-upper)
        df  :=    x   x   x   x   x
    */

    using namespace std::placeholders;

    auto st = std::bind(srhd::spherical_geometry_source_terms, _1, _2, M_PI / 2, gamma_law_index);
    auto xh = geometric::unit_vector_on_axis(1);
    auto bf = srhd::flux(initial_condition(run_config, front(state.vertices)), xh, gamma_law_index);
    auto dt = time_step(state);
    auto da = face_areas(state);
    auto dv = cell_volumes(state);
    auto dx = cell_spacing(state);
    auto x0 = cell_centers(state);
    auto p0 = state.conserved / dv | nd::map(recover_primitive()) | nd::to_shared();
    auto s0 = nd::zip(p0, x0) | nd::map(apply_to(st)) | nd::multiply(dv);
    auto gx = nd::zip(x0 | nd::adjacent_zip3(), p0 | nd::adjacent_zip3()) | nd::map(mara::plm_gradient(1.5)) | nd::extend_zeros() | nd::to_shared();
    auto pl = select(p0 + 0.5 * dx * gx, 0, 0, -1);
    auto pr = select(p0 - 0.5 * dx * gx, 0, 1);
    auto vf = nd::map((pl + pr) * 0.5, srhd::velocity_1) * 0.0 | nd::to_shared();
    auto pf = nd::zip(pl, pr);

    auto f0 = nd::zip(pf, vf)
    | nd::map(riemann_solver_for(xh))
    | nd::extend_uniform_lower(bf)
    | nd::extend_zero_gradient_upper()
    | nd::multiply(da)
    | nd::to_shared();

    auto df = f0 | nd::adjacent_diff() | nd::to_shared();
    auto q1 = state.conserved + (s0 - df) * dt | nd::to_shared();
    auto x1 = state.vertices  + (vf * dt | nd::extend_zero_gradient()) | nd::to_shared();

    return {
        state.iteration + 1,
        state.time + dt,
        x1,
        q1,
    };
}

auto advance(const mara::config_t& run_config)
{
    return [run_config] (state_t state)
    {
        return advance(run_config, state);
    };
}




//=============================================================================
void write_diagnostics(state_t state)
{
    std::printf("Write checkpoint %s\n", "sedov.h5");
    auto dv = cell_volumes(state);
    auto primitive = (state.conserved / dv) | nd::map(recover_primitive());
    auto file = h5::File("sedov.h5", "w");
    h5::write(file, "vertices", state.vertices | nd::adjacent_mean() | nd::to_shared());
    h5::write(file, "primitive", primitive | nd::to_shared());
}




//=============================================================================
auto should_continue(const mara::config_t& run_config)
{
    return [tfinal = run_config.get_double("tfinal")] (timed_state_pair_t p)
    {
        return control::last_state(p).time <= dimensional::unit_time(tfinal);
    };
}

void print_run_loop(timed_state_pair_t p)
{
    auto nz = double(size(control::this_state(p).vertices));
    auto us = control::microseconds_separating(p);
    auto s = control::this_state(p);
    std::printf("[%04lu] t=%6.5lf Mzps=%3.2lf\n", long(s.iteration), s.time.value, nz / us);
}

auto side_effects(mara::config_t& run_config)
{
    return [tfinal = run_config.get_double("tfinal")] (timed_state_pair_t p)
    {
        if (control::this_state(p).time < dimensional::unit_time(tfinal))
        {
            print_run_loop(p);            
        }

        if (control::this_state(p).time >= dimensional::unit_time(tfinal))
        {
            write_diagnostics(control::this_state(p));
        }
    };
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto run_config = config_template()
    .create()
    .update(mara::argv_to_string_map(argc, argv));

    auto simulation = seq::generate(initial_state(run_config), advance(run_config))
    | seq::pair_with(control::time_point_sequence())
    | seq::window()
    | seq::take_while(should_continue(run_config));

    mara::pretty_print(std::cout, "config", run_config);

    for (auto state_pair : simulation)
    {
        side_effects(run_config)(state_pair);
    }
    return 0;
}
