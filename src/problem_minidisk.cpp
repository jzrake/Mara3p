/**
 ==============================================================================
 Copyright 2020, Jonathan Zrake

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




#include <iostream>
#include "mara.hpp"
#include "parallel_mpi.hpp"




//=============================================================================
namespace mpi
{
    inline filtered_ostream_t& master_cout()
    {
        static auto os = filtered_ostream_t(std::cout, mpi::comm_world().rank() == 0);
        return os;
    }    
}




#if MARA_COMPILE_LOCALLY_ISOTHERMAL // <---------------------------------------
#include "app_config.hpp"
#include "app_control.hpp"
#include "app_filesystem.hpp"
#include "app_hdf5.hpp"
#include "app_problem.hpp"
#include "core_unit_test.hpp"
#include "core_util.hpp"
#include "mesh_amr.hpp"
#include "module_locally_isothermal.hpp"
#include "physics_two_body.hpp"




//=============================================================================
inline auto config_template()
{
    return mara::config_template()
    .item("block_size",            64, "number of zones per side")
    .item("depth",                  1, "number of levels in the mesh")
    .item("mesh_type",      "uniform", "mesh type: [uniform]")
    .item("tfinal",               0.1, "time to stop the simulation")
    .item("cfl",                  0.3, "CFL number [0.0, 1.0]")
    .item("plm",                  1.5, "PLM theta parameter [1.0, 2.0]")
    .item("cpi",                  0.1, "time between checkpoints")
    .item("rk",                     2, "Runge-Kutta order")
    .item("threads",                1, "number of threads to execute on")
    .item("fold",                   1, "number of encodings per time step batch")
    .item("restart",               "", "a checkpoint file to restart from")
    .item("outdir",               ".", "the directory where output files are written");
}




//=============================================================================
class locally_isothermal_problem_t : public mara::problem_base_t
{
public:
    using solution_t = modules::locally_isothermal::solution_t;

    dimensional::unit_time get_time_from_solution(std::any solution) const override
    {
        return std::any_cast<solution_t>(solution).time;
    }

    void write_solution(std::string filename, std::any solution) const override
    {
        h5::write(h5::File(filename, "r+"), "solution", std::any_cast<solution_t>(solution));
    }

    std::any read_solution(std::string filename) const override
    {
        auto solution = solution_t();
        h5::read(h5::File(filename, "r"), "solution", solution);
        return solution;
    }

    rational::number_t get_iteration_from_solution(std::any solution) const override
    {
        return std::any_cast<solution_t>(solution).iteration;
    }

    mpr::node_set_t get_computable_nodes(std::any solution) const override
    {
        auto sequence = std::any_cast<solution_t>(solution).conserved;
        auto state = start(sequence);
        auto nodes = mpr::node_set_t();

        while (state.has_value())
        {
            nodes.insert(obtain(sequence, state.value()).node());
            state = next(sequence, state.value());
        }
        return nodes;
    }
};




//=============================================================================
class minidisk_problem_t : public locally_isothermal_problem_t {
public:
    minidisk_problem_t(unsigned long zones) : zones(zones) {}
    unsigned long get_zone_count(std::any) const override { return zones; }
    std::string get_module_name() const override { return "minidisk"; }
private:
    unsigned long zones = 0;
};




//=============================================================================
// static auto create_mesh_geometry(const mara::config_t& cfg)
// {
//     return modules::locally_isothermal::mesh_geometry_t(1.0, cfg.get_int("block_size"));
// }

// static auto create_mesh_topology(const mara::config_t& cfg)
// {
//     auto mesh_type = cfg.get_string("mesh_type");
//     auto depth     = cfg.get_int("depth");

//     if (mesh_type == "uniform")
//     {
//         return bsp::uniform_mesh_tree<2>(depth);
//     }
//     throw std::runtime_error("create_mesh (invalid mesh type " + mesh_type + ")");
// }

// static auto initial_solution(const mara::config_t& cfg)
// {
//     auto softening_length = dimensional::unit_length(0.01);
//     auto omega_frame      = dimensional::unit_rate(0.0);

//     auto initial = [rs=softening_length, omega_frame] (numeric::array_t<dimensional::unit_length, 2> p)
//     {
//         auto [x, y] = as_tuple(p);
//         auto ph = two_body::potential(two_body::point_mass_t(), x, y, rs);
//         auto r0 = sqrt(x * x + y * y);
//         auto vp = sqrt(-ph) - omega_frame * r0;
//         auto vx = vp * (-y / r0);
//         auto vy = vp * ( x / r0);
//         return iso2d::primitive(1.0, vx, vy);
//     };

//     auto geom = create_mesh_geometry(cfg);
//     auto mesh = create_mesh_topology(cfg);
//     auto u    = modules::locally_isothermal::initial_conserved_tree(mesh, geom, initial);

//     return modules::locally_isothermal::solution_t{0, 0.0, u};
// }




//=============================================================================
void run_minidisk(int argc, const char* argv[])
{
    using namespace modules;
    using namespace std::placeholders;

    // auto temperature = [] (auto) { return dimensional::unit_specific_energy(1.0); };
    // auto viscosity   = modules::locally_isothermal::unit_viscosity(0.0);

    // auto args             = mara::argv_to_string_map(argc, argv);
    // auto cfg              = config_template().create().update(mara::restart_run_config(args)).update(args);
    // auto tfinal           = dimensional::unit_time(cfg.get_double("tfinal"));
    // auto rk_order         = cfg.get_int("rk");
    // auto plm_theta        = cfg.get_double("plm");
    // auto mesh             = create_mesh_topology(cfg);
    // auto geom             = create_mesh_geometry(cfg);
    // auto schedule         = mara::initial_schedule(cfg);
    // auto dx               = locally_isothermal::smallest_cell_size(mesh, geom);
    // auto zones            = locally_isothermal::total_zones(mesh, geom);
    // auto problem          = minidisk_problem_t(zones);
    // auto solution         = mara::initial_solution(problem, cfg, initial_solution(cfg));
    // auto vm               = dimensional::unit_velocity(5.0);
    // auto strategy         = mpr::mpi_multi_threaded_execution(cfg.get_int("threads"));
    // auto monitor          = mpr::execution_monitor_t();
    // auto dt               = dx / vm;
    // auto scheme           = std::bind(modules::locally_isothermal::updated_solution, _1, dt, geom, temperature, viscosity, plm_theta);

    // mpr::compute_all(solution.conserved, mpr::mpi_single_threaded_execution());
    // mara::pretty_print(mpi::master_cout(), "config", cfg);

    // problem.print_task_graph(cfg, schedule, control::advance_runge_kutta(scheme, rk_order, solution));        
    // problem.side_effects(cfg, schedule, solution, monitor);

    // while (solution.time < tfinal)
    // {
    //     solution = control::advance_runge_kutta(scheme, rk_order, solution);
    //     monitor = mpr::compute_all(solution.conserved, strategy);
    //     problem.side_effects(cfg, schedule, solution, monitor);
    // }
}

#else

void run_minidisk(int argc, const char* argv[])
{
    mpi::master_cout() << "[mara] program not available: recompile with MARA_COMPILE_LOCALLY_ISOTHERMAL = 1" << std::endl;
}

#endif // MARA_COMPILE_LOCALLY_ISOTHERMAL
