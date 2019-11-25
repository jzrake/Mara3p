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
#include "app_config.hpp"
#include "app_control.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_geometric.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_rational.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_sequence.hpp"
#include "physics_mhd.hpp"
#include "mesh_cartesian_3d.hpp"
#include "scheme_mhd.hpp"




//=============================================================================
static const unsigned zone_count = 64;
static const double gamma_law_index = 5. / 3;




//=============================================================================
struct solution_t
{
    rational::number_t                            iteration;
    dimensional::unit_time                        time;
    nd::shared_array<mhd::conserved_density_t, 3> conserved;
    nd::shared_array<mhd::unit_magnetic_field, 3> magnetic_flux_1;
    nd::shared_array<mhd::unit_magnetic_field, 3> magnetic_flux_2;
    nd::shared_array<mhd::unit_magnetic_field, 3> magnetic_flux_3;
};




//=============================================================================
using namespace std::placeholders;
using util::apply_to;
using util::compose;
using geometric::unit_vector_on;
using geometric::component;
using position_t            = geometric::euclidean_vector_t<dimensional::unit_length>;
using solution_with_tasks_t = std::pair<solution_t, control::task_t>;
using timed_state_pair_t    = control::timed_pair_t<solution_with_tasks_t>;




//=============================================================================
auto magnetic_flux(solution_t s) { return std::tuple(s.magnetic_flux_1, s.magnetic_flux_2, s.magnetic_flux_3); }
auto conserved(solution_t s) { return s.conserved; }
auto solution(solution_with_tasks_t p) { return p.first; }
auto tasks(solution_with_tasks_t p)   { return p.second; }




//=============================================================================
namespace h5 {

void write(const Group& group, const solution_t& solution)
{
    write(group, "iteration", solution.iteration);
    write(group, "time", solution.time);
    write(group, "conserved", solution.conserved);
    write(group, "magnetic_flux_1", solution.magnetic_flux_1);
    write(group, "magnetic_flux_2", solution.magnetic_flux_2);
    write(group, "magnetic_flux_3", solution.magnetic_flux_3);
}

}




//=============================================================================
mhd::primitive_t initial_primitive(position_t x, mhd::magnetic_field_vector_t b)
{
    auto x0 = geometric::euclidean_vector<dimensional::unit_length>(0.5, 0.5, 0.5);
    auto R = sqrt(length_squared(x - x0));
    auto d = R < dimensional::unit_length(0.125) ? 1.0 : 0.1;
    auto p = R < dimensional::unit_length(0.125) ? 1.0 : 0.125;
    return mhd::primitive(d, {}, p, b);
}

mhd::vector_potential_t initial_vector_potential(position_t x)
{
    return mhd::vector_potential_t{0.0, x.component_1().value * 4.0, 0.0};
}




//=============================================================================
solution_t initial_solution()
{
    auto N = zone_count;

    auto dl = pow<1>(dimensional::unit_length(1.0) / double(N));
    auto da = pow<2>(dimensional::unit_length(1.0) / double(N));

    auto A = [] (unsigned dir) { return compose(component(dir), initial_vector_potential); };
    auto p2c = std::bind(mhd::conserved_density, _1, gamma_law_index);

    auto xv = mesh::unit_lattice<dimensional::unit_length>(N + 1, N + 1, N + 1);
    auto xc = mesh::cell_positions(xv);
    auto [xe1, xe2, xe3] = mesh::edge_positions(xv);
    auto [ae1, ae2, ae3] = std::tuple(xe1 | nd::map(A(1)), xe2 | nd::map(A(2)), xe3 | nd::map(A(3)));
    auto [mf1, mf2, mf3] = mesh::solenoidal_difference(ae1 * dl, ae2 * dl, ae3 * dl);
    auto [bf1, bf2, bf3] = std::tuple(mf1 / da, mf2 / da, mf3 / da);
    auto bc = mesh::face_to_cell(bf1, bf2, bf3);
    auto pc = nd::zip(xc, bc) | nd::map(apply_to(initial_primitive));
    auto uc = pc | nd::map(p2c);

    return {
        rational::number(0),
        dimensional::unit_time(0.0),
        uc  | nd::to_shared(),
        bf1 | nd::to_shared(),
        bf2 | nd::to_shared(),
        bf3 | nd::to_shared(),
    };
}




//=============================================================================
solution_t advance(solution_t solution)
{
    auto [bf1, bf2, bf3] = magnetic_flux(solution);
    auto uc = conserved(solution);
    auto dl = dimensional::unit_length(1.0) / double(zone_count);
    auto up = advance(solution.time, uc, bf1, bf2, bf3, dl);

    return {
        solution.iteration + 1,
        std::get<0>(up),
        std::get<1>(up),
        std::get<2>(up),
        std::get<3>(up),
        std::get<4>(up),
    };
}

solution_with_tasks_t advance_app_state(solution_with_tasks_t state)
{
    auto cpi = dimensional::unit_time(0.01);

    return std::pair(
        advance(solution(state)),
        jump(tasks(state), solution(state).time, cpi));
}

solution_with_tasks_t initial_app_state(const mara::config_t& run_config)
{
    return std::pair(initial_solution(), control::task("write_checkpoint", 0.0));
}

bool should_continue(timed_state_pair_t state_pair)
{
    return solution(control::this_state(state_pair)).time < dimensional::unit_time(1.0);
}




//=============================================================================
void write_checkpoint(solution_t solution, long count)
{
    auto fname = util::format("mhd3d.%04lu.h5", count);
    std::printf("Write checkpoint %s\n", fname.data());
    write(h5::File(fname, "w"), solution);
}

void print_run_loop(timed_state_pair_t p)
{
    auto soln = solution(control::this_state(p));
    auto nz = size(soln.conserved);
    auto us = control::microseconds_separating(p);

    std::printf("[%07lu] t=%.6lf dt=%.2e Mzps=%.2lf\n",
        long(soln.iteration),
        soln.time.value,
        0.2,
        nz / us);
}

void side_effects(timed_state_pair_t p)
{
    auto this_soln = solution(control::this_state(p));
    auto last_soln = solution(control::last_state(p));
    auto this_task = tasks(control::this_state(p));
    auto last_task = tasks(control::last_state(p));

    if (this_task.count != last_task.count)
        write_checkpoint(last_soln, last_task.count);

    print_run_loop(p);
}




//=============================================================================
int main()
{
    auto run_config = mara::config_t();

    auto simulation = seq::generate(initial_app_state(run_config), advance_app_state)
    | seq::pair_with(control::time_point_sequence())
    | seq::window()
    | seq::take_while(should_continue);

    for (auto state_pair : simulation)
    {
        side_effects(state_pair);
    }
    return 0;
}
