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




//=============================================================================
static const auto gamma_law_index = 5. / 3;
static const auto num_cells = 8192;
static const auto final_time = 0.22;




//=============================================================================
struct state_t
{
    dimensional::unit_time time = 0.0;
    rational::number_t iteration = 0;
    nd::shared_array<euler::conserved_density_t, 1> conserved;
    std::chrono::time_point<std::chrono::high_resolution_clock> time_point;
};

using timed_state_pair_t = control::timed_pair_t<state_t>;




//=============================================================================
auto vertices(nd::uint cell_count)
{
    return nd::linspace(0.0, 1.0, cell_count + 1) | nd::construct<dimensional::unit_length>();
}

auto cell_centers(nd::uint cell_count)
{
    return vertices(cell_count) | nd::adjacent_mean();
}

auto cell_spacings(nd::uint cell_count)
{
    return vertices(cell_count) | nd::adjacent_diff();
}

euler::primitive_t initial_condition(dimensional::unit_length x)
{
    return x < dimensional::unit_length(0.5) ? euler::primitive(1.0, 1.0) : euler::primitive(0.125, 0.1);
}

auto riemann_solver_for(geometric::unit_vector_t nhat)
// -> std::function<euler::flux_vector_t(std::tuple<euler::primitive_t, euler::primitive_t>)>
{
    return [nhat] (auto t)
    {
        return euler::riemann_hlle(std::get<0>(t), std::get<1>(t), nhat, gamma_law_index);
    };
}

auto recover_primitive()
// -> std::function<euler::primitive_t(euler::conserved_density_t)>
{
    return [] (auto u) { return euler::recover_primitive(u, gamma_law_index); };
}




//=============================================================================
state_t initial_state()
{
    auto p0 = cell_centers(num_cells) | nd::map(initial_condition);
    auto u0 = p0 | nd::map([] (auto p) { return euler::conserved_density(p, gamma_law_index); });

    return {0.0, 0, u0 | nd::to_shared(), std::chrono::high_resolution_clock::now()};
}

state_t advance(state_t state)
{
    auto xh = geometric::unit_vector_on(1);
    auto dt = dimensional::unit_time(0.3 / num_cells);
    auto dx = cell_spacings(num_cells);
    auto u0 = state.conserved;
    auto p0 = u0 | nd::map(recover_primitive());
    auto f0 = p0 | nd::extend_zero_gradient() | nd::adjacent_zip() | nd::map(riemann_solver_for(xh)) | nd::to_shared();
    auto df = f0 | nd::adjacent_diff();
    auto u1 = u0 - df * dt / dx;

    return {
        state.time + dt,
        state.iteration + 1,
        u1 | nd::to_shared(),
        std::chrono::high_resolution_clock::now(),
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
        auto primitive = (control::this_state(p).conserved) | nd::map(recover_primitive());
        auto file = h5::File("euler.h5", "w");
        h5::write(file, "cell_centers", cell_centers(num_cells) | nd::to_shared());
        h5::write(file, "primitive", primitive | nd::to_shared());
    }
}

auto time_point_sequence()
{
    using namespace std::chrono;
    return seq::generate(high_resolution_clock::now(), [] (auto) { return high_resolution_clock::now(); });
}




//=============================================================================
int main()
{
    auto simulation = seq::generate(initial_state(), advance)
    | seq::pair_with(time_point_sequence())
    | seq::window()
    | seq::take_while(should_continue);

    for (auto state_pair : simulation)
    {
        side_effects(state_pair);
    }
    return 0;
}
