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
#include "scheme_mhd_v1.hpp"




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("nc",                64)   // number of cells on each axis
    .item("tf",              10.0)   // time to stop the simulation
    .item("cpi",              0.1)   // checkpoint (output) interval
    .item("rk",                 2)   // Runge-Kutta order (1, 2, or 3)
    .item("cfl",             0.33)   // courant number
    .item("plm",              1.5)   // PLM parameter
    .item("b0",               1.0)   // field strength parameter
    .item("setup", "spherical-blast");
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

solution_t weighted_sum(solution_t s, solution_t t, rational::number_t b)
{
    return {
        s.iteration  *         b + t.iteration *       (1 - b),
        s.time       * double(b) + t.time      * double(1 - b),
        (s.conserved * double(b) + t.conserved * double(1 - b)) | nd::to_shared(),
        (s.magnetic_flux_1 * double(b) + t.magnetic_flux_1 * double(1 - b)) | nd::to_shared(),
        (s.magnetic_flux_2 * double(b) + t.magnetic_flux_2 * double(1 - b)) | nd::to_shared(),
        (s.magnetic_flux_3 * double(b) + t.magnetic_flux_3 * double(1 - b)) | nd::to_shared(),
    };
}




//=============================================================================
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

    auto p = mara::primitive_array(
        solution.conserved,
        std::array{
            solution.magnetic_flux_1,
            solution.magnetic_flux_2,
            solution.magnetic_flux_3});

    write(group, "primitive", p);
}

}




//=============================================================================
std::tuple<mara::primitive_function_t, mara::vector_potential_function_t> initial_data_functions(const mara::config_t& cfg)
{
    auto setup = cfg.get_string("setup");
    auto b0 = cfg.get_double("b0");

    if (setup == "spherical-blast")
    {
        auto ip = [] (position_t x, mhd::magnetic_field_vector_t b)
        {
            auto x0 = position_t{0.5, 0.5, 0.5};
            auto R = sqrt(length_squared(x - x0));
            auto d = R < dimensional::unit_length(0.125) ? 1.0 : 0.1;
            auto p = R < dimensional::unit_length(0.125) ? 1.0 : 0.125;
            return mhd::primitive(d, {}, p, b);
        };
        auto iv = [b0] (position_t x)
        {
            return b0 * mhd::vector_potential_t{0.0, x.component_1().value, 0.0};
        };
        return std::tuple(ip, iv);
    }

    if (setup == "abc")
    {
        auto ip = [] (position_t x, mhd::magnetic_field_vector_t b)
        {
            return mhd::primitive(1.0, {}, 1.0, b);
        };
        auto iv = [b0] (position_t p)
        {
            auto k = 2.0 * M_PI / dimensional::unit_length(1.0);
            auto [A, B, C] = std::tuple(1.0, 1.0, 1.0);
            auto [x, y, z] = as_tuple(p);
            auto ax = A * std::sin(k * z) + C * std::cos(k * y);
            auto ay = B * std::sin(k * x) + A * std::cos(k * z);
            auto az = C * std::sin(k * y) + B * std::cos(k * x);
            return b0 * mhd::vector_potential_t{ax, ay, az};
        };
        return std::tuple(ip, iv);
    }

    throw std::runtime_error("mhd3d::initial_data_functions (no setup named " + setup + ")");
}




//=============================================================================
auto advance_solution(const mara::config_t& cfg)
{
    auto cfl = cfg.get_double("cfl");
    auto plm = cfg.get_double("plm");
    auto dl  = dimensional::unit_length(1.0) / double(cfg.get_int("nc"));

    return [cfl, plm, dl] (solution_t solution) -> solution_t
    {
        auto [bf1, bf2, bf3] = magnetic_flux(solution);
        auto bc = mara::local_periodic_boundary_extension();
        auto uc = conserved(solution);
        auto [t_, uc_, bf1_, bf2_, bf3_] = mara::advance_mhd_ct(solution.time, uc, bf1, bf2, bf3, dl, cfl, plm, *bc);
        return {solution.iteration + 1, t_, uc_, bf1_, bf2_, bf3_};
    };
}

auto advance_app_state(const mara::config_t& cfg)
{
    auto cpi          = dimensional::unit_time(1.0) * cfg.get_double("cpi");
    auto rk           = cfg.get_int("rk");
    auto advance_rk   = control::advance_runge_kutta(advance_solution(cfg), rk);

    return [cpi, advance_rk] (solution_with_tasks_t state)
    {
        auto next_soln = advance_rk(solution(state));

        return std::pair(
            next_soln,
            jump(tasks(state), next_soln.time, cpi));
    };
}

auto should_continue(const mara::config_t& cfg)
{
    auto tf = dimensional::unit_time(cfg.get_double("tf"));

    return [tf] (timed_state_pair_t state_pair)
    {
        return solution(control::last_state(state_pair)).time <= tf;
    };
}




//=============================================================================
solution_t initial_solution(const mara::config_t& cfg)
{
    auto nc = cfg.get_int("nc");
    auto [ip, iv] = initial_data_functions(cfg);
    auto [uc, bf1, bf2, bf3] = mara::construct_conserved(ip, iv, nc);
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
