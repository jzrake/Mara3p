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
#include "app_control.hpp"
#include "app_filesystem.hpp"
#include "app_hdf5_config.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"
#include "parallel_computable.hpp"
#include "parallel_computable_tree.hpp"
#include "parallel_thread_pool.hpp"
#include "physics_two_body.hpp"
#include "minidisk.hpp"
#include "minidisk_io.hpp"




/*
 * Todo list:
 * 
 * [x] eccentric orbits and non-equal mass binaries
 * [x] configurable frame rotation rate
 * [x] configurable Mach number
 * [x] compute (rather than estimate) allowed time step
 * [x] inclusion of viscous stress
 * [ ] variable size blocks in FMR mesh
 *     ( ) generalize extend method with prolongation
 *     ( ) flux correction
 *     ( ) upsampling in plotting script
 * [x] proper stepping of side effects
 * [x] checkpoint read / restart
 * [x] break up code into separate source files
 * [ ] VTK output option
 * [ ] time series of forces, work, orbital element evolution, etc.
 * [ ] computable: execution statistics to discover bottlenecks
 * [ ] computable: MPI execution strategy
 */




//=============================================================================
using namespace dimensional;
using namespace minidisk;
using namespace std::placeholders;
static const auto binary_period = unit_time(2 * M_PI);




//=============================================================================
template<typename ArrayType>
auto extend_x(bsp::shared_tree<pr::computable<ArrayType>, 4> __pc__, bsp::tree_index_t<2> block, bool dummy=false)
{
    if (dummy)
    {
        return map(value_at(__pc__, block), [] (auto x) { return x; }).immediate(true);
    }

    auto _pl_ = value_at(__pc__, prev_on(block, 0));
    auto _pc_ = value_at(__pc__, block);
    auto _pr_ = value_at(__pc__, next_on(block, 0));

    return pr::zip(_pl_, _pc_, _pr_) | pr::mapv([] (auto pl, auto pc, auto pr)
    {
        auto nx = shape(pc, 0);
        auto ny = shape(pc, 1);

        return nd::make_array(nd::indexing([pl, pc, pr, nx] (auto i, auto j)
        {
            if (i == 0)      return pl(nx - 1, j);
            if (i == nx + 1) return pr(0, j);
            return                  pc(i - 1, j);
        }), nd::uivec(nx + 2, ny)) | nd::to_shared();
    });
}

template<typename ArrayType>
auto extend_y(bsp::shared_tree<pr::computable<ArrayType>, 4> __pc__, bsp::tree_index_t<2> block, bool dummy=false)
{
    if (dummy)
    {
        return map(value_at(__pc__, block), [] (auto x) { return x; }).immediate(true);
    }

    auto _pl_ = value_at(__pc__, prev_on(block, 1));
    auto _pc_ = value_at(__pc__, block);
    auto _pr_ = value_at(__pc__, next_on(block, 1));

    return pr::zip(_pl_, _pc_, _pr_) | pr::mapv([] (auto pl, auto pc, auto pr)
    {
        auto nx = shape(pc, 0);
        auto ny = shape(pc, 1);

        return nd::make_array(nd::indexing([pl, pc, pr, ny] (auto i, auto j)
        {
            if (j == 0)      return pl(i, ny - 1);
            if (j == ny + 1) return pr(i, 0);
            return                  pc(i, j - 1);
        }), nd::uivec(nx, ny + 2)) | nd::to_shared();
    });
}




//=============================================================================
static auto initial_solution(const mara::config_t& cfg)
{
    auto solver_data = make_solver_data(cfg);
    auto restart = cfg.get_string("restart");

    if (! restart.empty())
    {
        auto result = solution_t{};
        auto h5f = h5::File(restart, "r");
        read(h5f, "solution", result);
        return result;
    }

    return solution_t{
        0,
        0.0,
        bsp::uniform_quadtree(solver_data.depth) | bsp::maps(std::bind(initial_conserved_array, solver_data, _1)),
    };
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
static auto smallest_cell_crossing_time(conserved_array_t uc, bsp::tree_index_t<2> block, solver_data_t solver_data)
{
    auto vm = nd::max(uc | nd::map(iso2d::recover_primitive) | nd::map(std::bind(iso2d::fastest_wavespeed, _1, unit_specific_energy(0.0))));
    auto dl = cell_size(block, solver_data);
    return dl / vm;
}

static auto smallest_cell_crossing_time(conserved_tree_t uc, solver_data_t solver_data)
{
    auto dt = unit_time(std::numeric_limits<double>::max());

    sink(indexify(uc), util::apply_to([&dt, solver_data] (auto block, auto uc)
    {
        dt = std::min(dt, smallest_cell_crossing_time(uc, block, solver_data));
    }));
    return dt;
}




//=============================================================================
static auto encode_step(
    rational::number_t iteration,
    unit_time time,
    bsp::shared_tree<pr::computable<conserved_array_t>, 4> conserved,
    unit_time dt,
    solver_data_t solver_data)
{
    auto skip_trans = solver_data.kinematic_viscosity == unit_viscosity(0.0);
    auto mesh   = indexes(conserved);
    auto _uc_   = conserved;
    auto _pc_   = _uc_ | bsp::maps(pr::map(minidisk::recover_primitive_array, "P"));
    auto _p_ext_x_ = mesh | bsp::maps([_pc_] (auto index) { return extend_x(_pc_, index).name("P-x"); });
    auto _p_ext_y_ = mesh | bsp::maps([_pc_] (auto index) { return extend_y(_pc_, index).name("P-y"); });
    auto _gradient_x_ = _p_ext_x_ | bsp::maps(pr::map(std::bind(estimate_gradient, _1, 0, 2.0), "Gx"));
    auto _gradient_y_ = _p_ext_y_ | bsp::maps(pr::map(std::bind(estimate_gradient, _1, 1, 2.0), "Gy"));

    auto _gradient_x_ext_x_ = mesh | bsp::map([g=_gradient_x_            ] (auto index) { return extend_x(g, index).name("Gx-x"); });
    auto _gradient_x_ext_y_ = mesh | bsp::map([g=_gradient_x_, skip_trans] (auto index) { return extend_y(g, index, skip_trans).name("Gx-y"); });
    auto _gradient_y_ext_x_ = mesh | bsp::map([g=_gradient_y_, skip_trans] (auto index) { return extend_x(g, index, skip_trans).name("Gy-x"); });
    auto _gradient_y_ext_y_ = mesh | bsp::map([g=_gradient_y_            ] (auto index) { return extend_y(g, index).name("Gy-y"); });

    auto _godunov_fluxes_x_ = zip(_p_ext_x_, _gradient_x_ext_x_, _gradient_y_ext_x_, mesh)
    | bsp::mapvs([solver_data] (auto pc, auto gl, auto gt, auto block)
    {
        return pr::zip(pc, gl, gt)
        | pr::mapv(std::bind(godunov_and_viscous_fluxes, _1, _2, _3, solver_data, block, 0), "Fx");
    });

    auto _godunov_fluxes_y_ = zip(_p_ext_y_, _gradient_y_ext_y_, _gradient_x_ext_y_, mesh)
    | bsp::mapvs([solver_data] (auto pc, auto gl, auto gt, auto block)
    {
        return pr::zip(pc, gl, gt)
        | pr::mapv(std::bind(godunov_and_viscous_fluxes, _1, _2, _3, solver_data, block, 1), "Fy");
    });

    auto u1 = zip(_uc_, _pc_, _godunov_fluxes_x_, _godunov_fluxes_y_, mesh)
    | bsp::mapvs([t=time, dt, solver_data] (auto uc, auto pc, auto fx, auto fy, auto block)
    {
        return pr::zip(uc, pc, fx, fy)
        | pr::mapv(std::bind(updated_conserved, _1, _2, _3, _4, t, dt, block, solver_data), "U");
    });

    return std::tuple(iteration + 1, time + dt, u1);
}




//=============================================================================
static auto step(solution_t solution, solver_data_t solver_data, int nfold, mara::ThreadPool& pool, bool print_only)
{
    auto [iter, time, uc] = std::tuple(
        solution.iteration,
        solution.time,
        solution.conserved | bsp::maps(pr::just<conserved_array_t>));

    auto cfl = solver_data.cfl * std::min(1.0, 0.01 + solution.time / unit_time(0.05));
    auto dt  = smallest_cell_crossing_time(solution.conserved, solver_data) * cfl;
    auto dt_est = cfl * smallest_cell_size(solver_data) / sqrt(two_body::G * unit_mass(0.5) / solver_data.softening_length);

    for (int i = 0; i < nfold; ++i)
    {
        std::tie(iter, time, uc) = encode_step(iter, time, uc, dt, solver_data);
    }

    auto master = computable(uc).name("_U_");

    if (print_only)
    {
        auto outfile = std::fopen("minidisk.dot", "w");
        pr::print_graph(outfile, master);
        std::fclose(outfile);
        return std::pair(solution, unit_scalar(0.0));
    }
    else
    {
        compute(master, pool.scheduler());
        return std::pair(solution_t{iter, time, master.value()}, dt / dt_est);
    }
}




//=============================================================================
static void side_effects(const mara::config_t& cfg, solution_t solution, schedule_t& schedule)
{
    if (solution.time >= schedule.checkpoint.next_due)
    {
        auto outdir = cfg.get_string("outdir");
        auto fname = util::format("%s/chkpt.%04d.h5", outdir.data(), schedule.checkpoint.count);

        mara::filesystem::require_dir(outdir);

        schedule.checkpoint.next_due += cfg.get_double("cpi") * binary_period;
        schedule.checkpoint.count += 1;

        auto h5f = h5::File(fname, "w");
        h5::write(h5f, "git_commit", std::string(GIT_COMMIT));
        h5::write(h5f, "run_config", cfg);
        h5::write(h5f, "solution", solution);
        h5::write(h5f, "schedule", schedule);
        std::printf("write %s\n", fname.data());
    }
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto args = mara::argv_to_string_map(argc, argv);
    auto cfg         = minidisk::config_template().create().update(restart_run_config(args)).update(args);
    auto solver_data = minidisk::make_solver_data(cfg);
    auto schedule    = initial_schedule(cfg);
    auto solution    = initial_solution(cfg);
    auto orbits      = cfg.get_double("orbits");
    auto fold        = cfg.get_int("fold");
    auto pool        = mara::ThreadPool(cfg.get_int("threads"));
    auto relative_dt = 0.0;

    step(solution, solver_data, fold, pool, true);
    mara::pretty_print(std::cout, "config", cfg);

    while (solution.time < orbits * binary_period)
    {
        side_effects(cfg, solution, schedule);

        auto [result, ticks] = control::invoke_timed(step, solution, solver_data, fold, pool, false);

        std::tie(solution, relative_dt) = result;

        auto ncells = std::pow(solver_data.block_size * (1 << solver_data.depth), 2);
        auto kzps = 1e6 * fold * ncells / std::chrono::duration_cast<std::chrono::nanoseconds>(ticks).count();

        std::printf("[%06ld] orbit=%lf dt(rel)=%0.2lf kzps=%lf\n", long(solution.iteration), solution.time.value / 2 / M_PI, relative_dt, kzps);
        std::fflush(stdout);
    }

    side_effects(cfg, solution, schedule);

    return 0;
}
