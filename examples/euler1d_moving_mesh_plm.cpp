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




#include <cstdio>
#include <chrono>
#include <tuple>
#include "app_control.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_sequence.hpp"
#include "physics_euler.hpp"
#include "scheme_plm_gradient.hpp"




//=============================================================================
static const auto gamma_law_index = 5. / 3;
static const auto xhat = geometric::unit_vector_t{{1.0, 0.0, 0.0}};
static const auto num_cells = 512;
static const auto final_time = 0.125;




//=============================================================================
struct state_t
{
    dimensional::unit_time time = 0.0;
    rational::number_t iteration = 0;
    nd::shared_array<dimensional::unit_length, 1> vertices;
    nd::shared_array<euler::conserved_t, 1> conserved;
};

using timed_state_pair_t = control::timed_pair_t<state_t>;




//=============================================================================
auto initial_vertices(nd::uint cell_count)
{
    return nd::linspace(0.0, 1.0, cell_count + 1) | nd::construct<dimensional::unit_length>();
}

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
    return (state.vertices | nd::adjacent_diff()) * dimensional::unit_area(1.0);
}

euler::primitive_t initial_condition(dimensional::unit_length x)
{
    // return x < dimensional::unit_length(0.5) ? euler::primitive(1.0, 1.0) : euler::primitive(0.125, 0.1);
    auto x0 = dimensional::unit_length(0.5);
    auto dx = dimensional::unit_length(0.05);
    auto d = 1.0 + std::exp(-std::pow((x - x0) / dx, 2.0));
    return euler::primitive(d, 1.0, 0.0, 0.0, 1.0);
}

auto riemann_solver_for(geometric::unit_vector_t nhat)
{
    return [nhat] (auto t)
    {
        auto [pf, vf] = t;
        auto [pl, pr] = pf;
        return euler::riemann_hlle(pl, pr, nhat, vf, gamma_law_index);
    };
}

auto recover_primitive()
{
    return [] (auto u) { return euler::recover_primitive(u, gamma_law_index); };
}




//=============================================================================
state_t initial_state()
{
    auto xv = initial_vertices(num_cells);
    auto dv = (xv | nd::adjacent_diff()) * dimensional::unit_area(1.0);
    auto p0 = xv | nd::adjacent_mean() | nd::map(initial_condition);
    auto u0 = p0 | nd::map([] (auto p) { return euler::conserved_density(p, gamma_law_index); });
    return {0.0, 0, xv | nd::to_shared(), (u0 * dv) | nd::to_shared()};
}

state_t advance(state_t state)
{
    auto dt = dimensional::unit_time(0.1 / num_cells);
    auto da = dimensional::unit_area(1.0);
    auto dv = cell_volumes(state);
    auto dx = cell_spacing(state);
    auto x0 = cell_centers(state);
    auto u0 = state.conserved / dv;
    auto p0 = u0 | nd::map(recover_primitive());
    auto gx = nd::zip(x0 | nd::adjacent_zip3(), p0 | nd::adjacent_zip3()) | nd::map(mara::plm_gradient(1.5)) | nd::extend_zero_gradient() | nd::to_shared();
    auto pl = select(p0 + 0.5 * dx * gx, 0, 0, -1);
    auto pr = select(p0 - 0.5 * dx * gx, 0, 1);
    auto pf = nd::zip(pl, pr);
    auto vf = 0.5 * (map(pl, euler::velocity_1) + map(pr, euler::velocity_1));

    auto f0 = zip(pf, vf) | nd::map(riemann_solver_for(xhat)) | nd::to_shared();
    auto df = f0 | nd::adjacent_diff() | nd::extend_zero_gradient();
    auto q1 = state.conserved - (df * dt * da) | nd::to_shared();
    auto x1 = state.vertices  + (vf * dt | nd::extend_zero_gradient()) | nd::to_shared();

    return {
        state.time + dt,
        state.iteration + rational::number(1),
        x1,
        q1,
    };
}

bool should_continue(timed_state_pair_t p)
{
    return control::last_state(p).time < dimensional::unit_time(final_time);
}

void print_run_loop(timed_state_pair_t p)
{
    auto us = control::microseconds_separating(p);
    auto s = control::this_state(p);
    std::printf("[%04lu] t=%6.5lf Mzps=%3.2lf\n", long(s.iteration), s.time.value, double(num_cells) / us);    
}

void side_effects(timed_state_pair_t p)
{
    print_run_loop(p);

    if (control::this_state(p).time >= dimensional::unit_time(final_time))
    {
        auto dv = cell_volumes(control::this_state(p));
        auto primitive = (control::this_state(p).conserved / dv) | nd::map(recover_primitive());
        auto file = h5::File("euler.h5", "w");
        h5::write(file, "cell_centers", control::this_state(p).vertices | nd::adjacent_mean() | nd::to_shared());
        h5::write(file, "primitive", primitive | nd::to_shared());
    }
}




//=============================================================================
int main()
{
    auto simulation = seq::generate(initial_state(), advance)
    | seq::pair_with(control::time_point_sequence())
    | seq::window()
    | seq::take_while(should_continue);

    for (auto state_pair : simulation)
    {
        side_effects(state_pair);
    }
    return 0;
}
