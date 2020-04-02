/**
 ==============================================================================
 Copyright 2019 - 2020, Jonathan Zrake

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
#include "app_filesystem.hpp"
#include "app_hdf5.hpp"
#include "app_problem.hpp"
#include "core_unit_test.hpp"
#include "mara.hpp"
#include "mesh_amr.hpp"
#include "module_euler.hpp"
#include "parallel_mpi.hpp"




/*
 * Todo list:
 * 
 * [x] put in runtime parameters
 * [x] proper side effects
 * [x] write ensure_valid_quadtree function
 * [x] refactor problem / module / scheme code to reduce boilerplate
 * [x] improve FMR scaling (cache computed blocks maybe)
 * [ ] check breaking of conservation laws at refinement boundaries
 * [ ] write the flux correction
 * 
 */




//=============================================================================
inline auto config_template()
{
    return mara::config_template()
    .item("block_size",            64, "number of zones per side")
    .item("depth",                  1, "number of levels in the mesh")
    .item("mesh_type",      "uniform", "mesh type: [uniform|nested]")
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
static auto create_mesh_geometry(const mara::config_t& cfg)
{
    return modules::euler2d::mesh_geometry_t(1.0, cfg.get_int("block_size"));
}

static auto create_mesh_topology(const mara::config_t& cfg)
{
    auto geom      = create_mesh_geometry(cfg);
    auto mesh_type = cfg.get_string("mesh_type");
    auto depth     = cfg.get_int("depth");

    if (mesh_type == "uniform")
    {
        return bsp::uniform_quadtree(depth);
    }
    if (mesh_type == "nested")
    {
        return amr::valid_quadtree(bsp::quadtree([geom] (auto block)
        {
            auto p = geom.block_centroid(block);
            auto r = sqrt(sum(p * p));
            return r < dimensional::unit_length(1.0) / double(block.level);
        }, depth));
    }
    throw std::runtime_error("create_mesh (invalid mesh type " + mesh_type + ")");
}

static auto initial_solution(const mara::config_t& cfg)
{
    auto shocktube_2d = [] (auto x)
    {
        if (sqrt(sum(x * x)) < dimensional::unit_length(0.2))
        {
            return euler::primitive(1.0, 0.0, 0.0, 0.0, 1.0);
        }
        else
        {
            return euler::primitive(0.1, 0.0, 0.0, 0.0, 0.125);
        }
    };

    auto gamma_law_index  = 5. / 3;
    auto geom             = create_mesh_geometry(cfg);
    auto mesh             = create_mesh_topology(cfg);
    auto u                = mpr::compute_all(modules::euler2d::initial_conserved_tree(mesh, geom, shocktube_2d, gamma_law_index), 1);

    return modules::euler2d::solution_t{0, 0.0, u};
}




//=============================================================================
class euler2d_demo_problem_t : public mara::problem_base_t
{
public:
    dimensional::unit_time get_time_from_solution(std::any solution) const override
    {
        return std::any_cast<modules::euler2d::solution_t>(solution).time;
    }

    void write_solution(std::string filename, std::any solution) const override
    {
        h5::write(h5::File(filename, "r+"), "solution", std::any_cast<modules::euler2d::solution_t>(solution));
    }

    std::any read_solution(std::string filename) const override
    {
        auto solution = modules::euler2d::solution_t();
        h5::read(h5::File(filename, "r+"), "solution", solution);
        return solution;
    }
};




//=============================================================================
static void run_euler2d(int argc, const char* argv[])
{
    using namespace modules::euler2d;
    using namespace std::placeholders;

    auto args             = mara::argv_to_string_map(argc, argv);
    auto cfg              = config_template().create().update(mara::restart_run_config(args)).update(args);
    auto problem          = euler2d_demo_problem_t();
    auto tfinal           = dimensional::unit_time(cfg.get_double("tfinal"));
    auto rk_order         = cfg.get_int("rk");
    auto plm_theta        = cfg.get_double("plm");
    auto gamma_law_index  = 5. / 3;
    auto mesh             = create_mesh_topology(cfg);
    auto geom             = create_mesh_geometry(cfg);
    auto schedule         = mara::initial_schedule(cfg);
    auto solution         = mara::initial_solution(problem, cfg, initial_solution(cfg));
    auto dx               = smallest_cell_size(mesh, geom);
    auto kz               = total_cells(mesh, geom) / 1e3;
    auto vm               = dimensional::unit_velocity(5.0);

    mara::pretty_print(std::cout, "config", cfg);

    while (solution.time < tfinal)
    {
        problem.side_effects(cfg, schedule, solution);

        auto dt = dx / vm;
        auto scheme = std::bind(updated_solution, _1, dt, geom, plm_theta, gamma_law_index);
        solution = control::advance_runge_kutta(scheme, rk_order, solution);

        auto start = std::chrono::high_resolution_clock::now();
        mpr::compute_all(solution.conserved, cfg.get_int("threads"));
        auto ticks = std::chrono::high_resolution_clock::now() - start;
        std::printf("[%06ld] t=%.4lf kzps=%.2lf\n", long(solution.iteration), solution.time.value, kz / (ticks.count() * 1e-9));
    }
    problem.side_effects(cfg, schedule, solution);
}




//=============================================================================
int test_app();
int test_core();
int test_mesh();
int test_model();
int test_physics();




//=============================================================================
static void run_all_tests(int argc, const char* argv[])
{
    start_unit_tests();
    test_app();
    test_core();
    test_mesh();
    test_model();
    test_physics();
    report_test_results();
}




//=============================================================================
int main(int argc, const char* argv[])
{
    if (argc == 1)
    {
        std::printf("usage: mara [command] <key-val parameters>\n");
        return 0;
    }

    auto mpi_session = mpi::Session();
    auto command = std::string(argv[1]);

    if (false) {}
    else if (command == "test")    run_all_tests(argc - 1, argv + 1);
    else if (command == "euler2d") run_euler2d(argc - 1, argv + 1);
    else std::printf("unknown command: %s\n", command.data());

    return 0;
}
