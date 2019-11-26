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
#include <iostream>
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
auto config_template()
{
    return mara::config_template()
    .item("block_size",            64)   // number of cells on each axis
    .item("tfinal",              10.0)   // time to stop the simulation
    .item("cpi",                  0.1)   // checkpoint (output) interval
    .item("rk_order",               2)   // Runge-Kutta order (1 or 2)
    .item("cfl",                  0.5)   // courant number
    .item("plm_theta",            1.5)   // PLM parameter
    .item("setup",  std::string("spherical-blast"));
}




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
std::tuple<mara::primitive_function_t, mara::vector_potential_function_t> initial_data_functions(const mara::config_t& cfg)
{
    auto setup = cfg.get_string("setup");

    if (setup == "spherical-blast")
    {
        auto ip = [] (position_t x, mhd::magnetic_field_vector_t b)
        {
            auto x0 = geometric::euclidean_vector<dimensional::unit_length>(0.5, 0.5, 0.5);
            auto R = sqrt(length_squared(x - x0));
            auto d = R < dimensional::unit_length(0.125) ? 1.0 : 0.1;
            auto p = R < dimensional::unit_length(0.125) ? 1.0 : 0.125;
            return mhd::primitive(d, {}, p, b);
        };
        auto iv = [] (position_t x)
        {
            return mhd::vector_potential_t{0.0, x.component_1().value * 4.0, 0.0};
        };
        return std::tuple(ip, iv);
    }
    throw std::runtime_error("mhd3d::initial_data_functions (no setup named " + setup + ")");
}




//=============================================================================
auto advance_solution(const mara::config_t& cfg)
{
    auto cfl = cfg.get_double("cfl");

    return [cfl] (solution_t solution, dimensional::unit_length dl) -> solution_t
    {
        auto [bf1, bf2, bf3] = magnetic_flux(solution);
        auto bc = mara::local_periodic_boundary_extension();
        auto uc = conserved(solution);
        auto [t_, uc_, bf1_, bf2_, bf3_] = mara::advance_mhd_ct(solution.time, uc, bf1, bf2, bf3, dl, cfl, *bc);
        return {solution.iteration + 1, t_, uc_, bf1_, bf2_, bf3_};
    };
}

auto advance_app_state(const mara::config_t& cfg)
{
    auto cpi          = dimensional::unit_time  (1.0) * cfg.get_double("cpi");
    auto cell_size    = dimensional::unit_length(1.0) / double(cfg.get_int("block_size"));
    auto advance_soln = advance_solution(cfg);

    return [cpi, cell_size, advance_soln] (solution_with_tasks_t state)
    {
        return std::pair(
            advance_soln(solution(state), cell_size),
            jump(tasks(state), solution(state).time, cpi));
    };
}

auto should_continue(const mara::config_t& cfg)
{
    auto final_time = dimensional::unit_time(cfg.get_double("tfinal"));

    return [final_time] (timed_state_pair_t state_pair)
    {
        return solution(control::this_state(state_pair)).time < final_time;
    };
}




//=============================================================================
solution_t initial_solution(const mara::config_t& cfg)
{
    auto block_size = cfg.get_int("block_size");
    auto [ip, iv] = initial_data_functions(cfg);
    auto [uc, bf1, bf2, bf3] = mara::construct_conserved(ip, iv, block_size);
    return {rational::number(0), dimensional::unit_time(0.0), uc, bf1, bf2, bf3};
}

solution_with_tasks_t initial_app_state(const mara::config_t& cfg)
{
    return std::pair(initial_solution(cfg), control::task("write_checkpoint", 0.0));
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
    auto dt = soln.time - solution(control::last_state(p)).time;

    std::printf("[%07lu] t=%.6lf dt=%.2e Mzps=%.2lf\n",
        long(soln.iteration),
        soln.time.value,
        dt.value,
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
int main(int argc, const char* argv[])
{
    auto cfg = config_template()
    .create()
    .update(mara::argv_to_string_map(argc, argv));

    mara::pretty_print(std::cout, "config", cfg);

    auto simulation = seq::generate(initial_app_state(cfg), advance_app_state(cfg))
    | seq::pair_with(control::time_point_sequence())
    | seq::window()
    | seq::take_while(should_continue(cfg));

    for (auto state_pair : simulation)
    {
        side_effects(state_pair);
    }
    return 0;
}
