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




#include <fstream>
#include <iostream>
#include "mara.hpp"
#include "app_control.hpp"
#include "app_filesystem.hpp"
#include "app_hdf5_config.hpp"
#include "app_serial_ndarray.hpp"
#include "app_serial_numeric_tuple.hpp"
#include "app_serial_std_tuple.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"
#include "minidisk.hpp"
#include "minidisk_io.hpp"
#include "parallel_computable.hpp"
#include "parallel_mpi.hpp"
#include "physics_two_body.hpp"
#ifndef MARA_GIT_COMMIT
#define MARA_GIT_COMMIT ""
#endif




/*
 * Todo list:
 * 
 * [x] eccentric orbits and non-equal mass binaries
 * [x] configurable frame rotation rate
 * [x] configurable Mach number
 * [ ] compute (rather than estimate) allowed time step
 * [x] inclusion of viscous stress
 * [ ] variable size blocks in FMR mesh
 *     ( ) generalize extend method with prolongation
 *     ( ) flux correction
 *     ( ) upsampling in plotting script
 * [ ] RK time stepping
 * [x] proper stepping of side effects
 * [x] checkpoint read / restart
 * [x] break up code into separate source files
 * [ ] VTK output option
 * [ ] time series of forces, work, orbital element evolution, etc.
 * [ ] computable: execution statistics to discover bottlenecks
 * [x] computable: MPI execution strategy
 */




//=============================================================================
using namespace dimensional;
using namespace minidisk;
using namespace std::placeholders;
static const auto binary_period = unit_time(2 * M_PI);




//=============================================================================
template<typename T, std::size_t Ratio>
mpr::node_set_t computable_tree_nodes(bsp::shared_tree<T, Ratio> tree)
{
    auto nodes = mpr::node_set_t();
    sink(tree, [&nodes] (auto c) { nodes.insert(c.node()); });
    return nodes;
}




//=============================================================================
template<typename ArrayType>
auto extend(bsp::shared_tree<mpr::computable<ArrayType>, 4> tree, bsp::tree_index_t<2> block)
{
    auto c11 = value_at(tree, block);
    auto c00 = value_at(tree, prev_on(prev_on(block, 0), 1));
    auto c02 = value_at(tree, prev_on(next_on(block, 0), 1));
    auto c20 = value_at(tree, next_on(prev_on(block, 0), 1));
    auto c22 = value_at(tree, next_on(next_on(block, 0), 1));
    auto c01 = value_at(tree, prev_on(block, 0));
    auto c21 = value_at(tree, next_on(block, 0));
    auto c10 = value_at(tree, prev_on(block, 1));
    auto c12 = value_at(tree, next_on(block, 1));

    return mpr::zip(c00, c01, c02, c10, c11, c12, c20, c21, c22)
    | mpr::mapv([] (auto c00, auto c01, auto c02, auto c10, auto c11, auto c12, auto c20, auto c21, auto c22)
    {
        auto nx = shape(c11, 0);
        auto ny = shape(c11, 1);
        auto cs = std::array{
            std::array{c00, c01, c02},
            std::array{c10, c11, c12},
            std::array{c20, c21, c22},
        };

        return nd::make_array(nd::indexing([cs, nx, ny] (auto i, auto j)
        {
            auto bi = i < 2 ? 0 : (i >= nx + 2 ? 2 : 1);
            auto bj = j < 2 ? 0 : (j >= ny + 2 ? 2 : 1);
            auto ii = i < 2 ? i - 2 + nx  : (i >= nx + 2 ? i - 2 - nx : i - 2);
            auto jj = j < 2 ? j - 2 + ny  : (j >= ny + 2 ? j - 2 - ny : j - 2);
            return cs[bi][bj](ii, jj);
        }), nd::uivec(nx + 4, ny + 4)) | nd::to_shared();
    });
}




//=============================================================================
static auto initial_solution(const mara::config_t& cfg)
{
    if (auto restart = cfg.get_string("restart"); ! restart.empty())
    {
        auto result = solution_t{};
        auto h5f = h5::File(restart, "r");
        read(h5f, "solution", result);
        return result;
    }
    auto sd = make_solver_data(cfg);
    auto mesh = bsp::uniform_quadtree(sd.depth);

    return solution_t{
        0,
        0.0,
        mesh | bsp::maps([sd] (auto b) { return mpr::from(std::bind(initial_conserved_array, sd, b)).name("U"); })};
}

static auto initial_schedule(const mara::config_t& cfg)
{
    auto restart = cfg.get_string("restart");

    if (! restart.empty())
    {
        auto result = minidisk::schedule_t{};
        auto h5f = h5::File(restart, "r");
        read(h5f, "schedule", result);
        return result;
    }
    return minidisk::schedule_t();
}

static auto restart_run_config(const mara::config_string_map_t& args)
{
    if (args.count("restart"))
    {
        auto file = h5::File(args.at("restart"), "r");
        return h5::read<mara::config_parameter_map_t>(file, "run_config");
    }
    return mara::config_parameter_map_t{};
}




//=============================================================================
static auto encode_substep(solution_t solution, unit_time dt, solver_data_t solver_data)
{
    auto updated_u = [solution, dt, solver_data] (auto b)
    {
        return [b, t=solution.time, dt, solver_data] (auto uc, auto pe)
        {
            return minidisk::updated_conserved(uc, pe, t, dt, b, solver_data);
        };
    };

    auto mesh  = indexes(solution.conserved);
    auto uc    = solution.conserved;
    auto pc    = uc   | bsp::maps(mpr::map(minidisk::recover_primitive_array, "P"));
    auto pe    = mesh | bsp::maps([pc] (auto b) { return extend(pc, b); });
    auto u1    = mesh | bsp::maps([uc, pe, updated_u] (auto b) { return zip(value_at(uc, b), value_at(pe, b)) | mpr::mapv(updated_u(b)); });

    return solution_t{solution.iteration + 1, solution.time + dt, u1};
}




//=============================================================================
static auto compute(solution_t solution, unsigned threads=0)
{
    mpr::compute_mpi(computable_tree_nodes(solution.conserved), threads);
    return solution;
}

static auto encode_step(solution_t solution, solver_data_t solver_data, int nfold)
{
    auto vmax = sqrt(two_body::G * unit_mass(0.5) / solver_data.softening_length);
    auto cfl  = solver_data.cfl * std::min(1.0, 0.01 + solution.time / unit_time(0.05));
    auto dt   = smallest_cell_size(solver_data) / vmax * cfl;

    for (int i = 0; i < nfold; ++i)
    {
        solution = encode_substep(solution, dt, solver_data);
    }
    return solution;
}

static auto step(solution_t solution, solver_data_t solver_data, int nfold, unsigned threads=0)
{
    return compute(encode_step(solution, solver_data, nfold), threads);
}

static void print_graph(const mara::config_t& cfg, solution_t solution, solver_data_t solver_data, int nfold)
{
    auto outdir = cfg.get_string("outdir");
    auto fname  = util::format("%s/task_graph.dot", outdir.data());
    auto next   = encode_step(solution, solver_data, nfold);
    auto nodes  = computable_tree_nodes(next.conserved);

    if (mpi::comm_world().rank() == 0)
    {
        std::printf("write %s\n", fname.data());
        mara::filesystem::require_dir(outdir);
        auto outf = std::ofstream(fname);
        mpr::print_graph(outf, nodes);        
    }
}




//=============================================================================
static void side_effects(const mara::config_t& cfg, solution_t solution, schedule_t& schedule)
{
    if (solution.time >= schedule.checkpoint.next_due)
    {
        auto outdir = cfg.get_string("outdir");
        auto fname = util::format("%s/chkpt.%04d.h5", outdir.data(), schedule.checkpoint.count);

        schedule.checkpoint.next_due += cfg.get_double("cpi") * binary_period;
        schedule.checkpoint.count += 1;

        mpi::comm_world().invoke([&] ()
        {
            if (mpi::comm_world().rank() == 0)
            {
                mara::filesystem::require_dir(outdir);

                auto h5f = h5::File(fname, "w");
                h5::write(h5f, "git_commit", std::string(MARA_GIT_COMMIT));
                h5::write(h5f, "run_config", cfg);
                h5::write(h5f, "schedule", schedule);
                std::printf("write %s\n", fname.data());
            }
            h5::write(h5::File(fname, "r+"), "solution", solution);
        });
    }
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto session     = mpi::Session();
    auto args        = mara::argv_to_string_map(argc, argv);
    auto cfg         = minidisk::config_template().create().update(restart_run_config(args)).update(args);
    auto solver_data = minidisk::make_solver_data(cfg);
    auto schedule    = initial_schedule(cfg);
    auto solution    = compute(initial_solution(cfg), 1);
    auto orbits      = cfg.get_double("orbits");
    auto fold        = cfg.get_int("fold");
    auto threads     = cfg.get_int("threads");
    auto ticks       = std::chrono::high_resolution_clock::now().time_since_epoch();

    if (mpi::comm_world().rank() == 0)
    {
        mara::pretty_print(std::cout, "config", cfg);
        std::printf("\trun on %d MPI processes and %d threads (%d compute units)\n\n",
            mpi::comm_world().size(),
            threads,
            threads * mpi::comm_world().size());
    }
    print_graph(cfg, solution, solver_data, fold);        

    while (solution.time < orbits * binary_period)
    {
        side_effects(cfg, solution, schedule);
        std::tie(solution, ticks) = control::invoke_timed(step, solution, solver_data, fold, threads);

        auto ncells = std::pow(solver_data.block_size * (1 << solver_data.depth), 2);
        auto kzps = 1e6 * fold * ncells / ticks.count();

        if (mpi::comm_world().rank() == 0)
        {
            std::printf("[%06ld] orbit=%lf kzps=%lf\n", long(solution.iteration), solution.time.value / 2 / M_PI, kzps);
            std::fflush(stdout);
        }
    }

    side_effects(cfg, solution, schedule);
    return 0;
}
