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
#include "app_config.hpp"
#include "app_control.hpp"
#include "app_filesystem.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_config.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_rational.hpp"
#include "core_bqo_tree.hpp"
#include "core_bsp_tree.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_sequence.hpp"
#include "core_util.hpp"
#include "parallel_computable.hpp"
#include "parallel_thread_pool.hpp"
#include "physics_iso2d.hpp"
#include "physics_two_body.hpp"
#include "scheme_plm_gradient.hpp"




/*
 * Todo list:
 * 
 * [x] eccentric orbits and non-equal mass binaries
 * [ ] configurable frame rotation rate
 * [ ] configurable Mach number
 * [ ] inclusion of viscous stress
 * [ ] variable size blocks in FMR mesh
 *     ( ) generalize extend method with prolongation
 *     ( ) flux correction
 *     ( ) upsampling in plotting script
 * [ ] proper stepping of side effects
 * [ ] checkpoint read / restart
 * [ ] VTK output option
 * [ ] time series of forces, work, orbital element evolution, etc.
 * [ ] computable execution statistics to discover bottlenecks
 */




//=============================================================================
using namespace dimensional;
using namespace std::placeholders;

static const auto mach_number   = 10.0;
static const auto omega_frame   = unit_rate(1.0);
static const auto binary_period = unit_time(2 * M_PI);




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("block_size",            32)   // number of zones per side
    .item("depth",                  2)   // number of levels in the mesh
    .item("domain_radius",        1.5)   // half-size of the domain
    .item("eccentricity",         0.0)   // binary orbital eccentricity
    .item("mass_ratio",           1.0)   // binary mass ratio M2 / M1
    .item("softening_length",    0.01)   // gravitational softening length
    .item("sink_radius",         0.01)   // radius of mass sinks
    .item("sink_rate",            1e3)   // rate of mass and momentum removal at sinks
    .item("buffer_rate",          1e2)   // maximum rate of buffer driving
    .item("buffer_scale",         0.2)   // buffer onset distance
    .item("cfl",                  0.5)   // Courant number
    .item("orbits",             100.0)   // time to stop the simulation
    .item("cpi",                 1000)   // checkpoint interval
    .item("threads",                1)   // number of threads to execute on
    .item("fold",                   1)   // number of encodings per time step batch
    .item("restart", std::string(""))    // a checkpoint file to restart from
    .item("outdir",  std::string("."));  // the directory where output files are written
}

struct solver_data_t
{
    int                 block_size;
    int                 depth;
    unit_length         domain_radius;
    unit_length         softening_length;
    unit_scalar         eccentricity;
    unit_scalar         mass_ratio;
    unit_length         sink_radius;
    unit_rate           sink_rate;
    unit_rate           buffer_rate;
    unit_length         buffer_scale;
    unit_scalar         cfl;
};

auto make_solver_data(const mara::config_t& cfg)
{
    return solver_data_t{
        cfg.get_int("block_size"),
        cfg.get_int("depth"),
        cfg.get_double("domain_radius"),
        cfg.get_double("softening_length"),
        cfg.get_double("eccentricity"),
        cfg.get_double("mass_ratio"),
        cfg.get_double("sink_radius"),
        cfg.get_double("sink_rate"),
        cfg.get_double("buffer_rate"),
        cfg.get_double("buffer_scale"),
        cfg.get_double("cfl"),
    };
}




//=============================================================================
using conserved_array_t = nd::shared_array<iso2d::conserved_density_t, 2>;
using primitive_array_t = nd::shared_array<iso2d::primitive_t, 2>;
using godunov_f_array_t = nd::shared_array<iso2d::flux_vector_t, 2>;
using conserved_tree_t  = bsp::shared_tree<conserved_array_t, 4>;

struct solution_t
{
    rational::number_t      iteration;
    dimensional::unit_time  time;
    conserved_tree_t        conserved;
};




//=============================================================================
namespace h5 {

void write(const Group& group, std::string name, const solution_t& solution)
{
    auto sgroup = group.require_group(name);
    write(sgroup, "iteration", solution.iteration);
    write(sgroup, "time", solution.time);

    auto ugroup = sgroup.require_group("conserved");

    sink(indexify(solution.conserved), util::apply_to([&ugroup] (auto index, auto value)
    {
        write(ugroup, format_tree_index(index), value);
    }));
}

}




//=============================================================================
namespace bsp {

template<typename ValueType, uint Ratio>
auto computable(shared_tree<pr::computable<ValueType>, Ratio> tree)
{
    auto nodes = pr::computable_node_t::set_t();
    sink(tree, [&nodes] (auto& b) { nodes.insert(b.node()); });

    return pr::computable_t<shared_tree<ValueType, Ratio>>([tree] ()
    {
        return tree | bsp::maps([] (auto v) { return v.value(); });
    }, nodes);
}

}




//=============================================================================
template<typename KeyType, typename ValueType>
class memoizer_t
{
public:
    memoizer_t(std::function<ValueType(KeyType)> function) : function(function) {}

    ValueType operator()(const KeyType& arg_tuple)
    {
        auto lock = std::lock_guard<std::mutex>(mutex);

        if (! cache.count(arg_tuple))
        {
            cache[arg_tuple] = function(arg_tuple);
        }
        return cache[arg_tuple];
    }

    template<typename... Args>
    ValueType operator()(Args... args)
    {
        return this->operator()(std::tuple(args...));
    }

private:
    std::mutex mutex;
    std::function<ValueType(KeyType)> function;
    std::map<KeyType, ValueType> cache;
};

template<typename... Args, typename Function>
auto memoizer(Function function)
{
    using arg_type = std::tuple<Args...>;
    using res_type = std::invoke_result_t<Function, Args...>;

    return memoizer_t<arg_type, res_type>([function] (arg_type args)
    {
        return std::apply(function, args);
    });
}




//=============================================================================
auto initial_primitive(unit_length softening_length)
{
    return [rs=softening_length] (numeric::array_t<unit_length, 2> p)
    {
        auto [x, y] = as_tuple(p);
        auto ph = two_body::potential(two_body::point_mass_t(), x, y, rs);
        auto r0 = sqrt(x * x + y * y);
        auto vp = sqrt(-ph) - omega_frame * r0;
        auto vx = vp * (-y / r0);
        auto vy = vp * ( x / r0);
        return iso2d::primitive(1.0, vx, vy);
    };
}

auto sound_speed_squared(numeric::array_t<unit_length, 2> p, unit_length softening_length)
{
    auto [x, y] = as_tuple(p);
    auto ph = two_body::potential(two_body::point_mass_t(), x, y, softening_length);
    return -ph / mach_number / mach_number;
}

auto buffer_rate_field(unit_length domain_radius, unit_length buffer_scale, unit_rate buffer_rate)
{
    return [=] (numeric::array_t<unit_length, 2> p)
    {
        auto r = sqrt(sum(p * p));
        auto y = (r - domain_radius) / buffer_scale;
        return 0.5 * buffer_rate * (1.0 + std::tanh(y));
    };
}

auto sink_rate_field(solver_data_t solver_data, numeric::array_t<unit_length, 2> sink_position)
{
    return [solver_data, sink_position] (numeric::array_t<unit_length, 2> p)
    {
        auto r6 = pow<3>(sum((p - sink_position) * (p - sink_position)));
        auto s6 = pow<6>(solver_data.sink_radius);
        return solver_data.sink_rate * std::exp(-r6 / s6);
    };
}

auto coriolis_term(geometric::euclidean_vector_t<unit_velocity> v)
{
    auto [vx, vy, vz] = as_tuple(v);
    return 2.0 * omega_frame * numeric::array(vy, -vx);
}

auto centrifugal_term(numeric::array_t<unit_length, 2> p)
{
    auto [x, y] = as_tuple(p);
    return omega_frame * omega_frame * numeric::array(x, y);
}




//=============================================================================
auto vertex_positions(uint n, unit_length dr, bsp::tree_index_t<2> block)
{
    auto [i0, j0] = as_tuple(block.coordinates);
    auto [i1, j1] = std::tuple(i0 + 1, j0 + 1);
    auto dl = double(1 << block.level);

    auto x0 = dr * (-1 + 2 * i0 / dl);
    auto x1 = dr * (-1 + 2 * i1 / dl);
    auto y0 = dr * (-1 + 2 * j0 / dl);
    auto y1 = dr * (-1 + 2 * j1 / dl);

    auto xv = nd::linspace(x0, x1, n + 1);
    auto yv = nd::linspace(y0, y1, n + 1);

    return std::tuple(xv, yv);
}

auto face_coordinates(bsp::uint n, unit_length dr, bsp::tree_index_t<2> block, bsp::uint axis)
{
    auto [xv, yv] = vertex_positions(n, dr, block);
    auto xc = xv | nd::adjacent_mean(0);
    auto yc = yv | nd::adjacent_mean(0);
    auto vec2 = nd::map(util::apply_to([] (auto x, auto y) { return numeric::array(x, y); }));

    switch (axis)
    {
        case 0: return nd::cartesian_product(xv, yc) | vec2 | nd::to_shared();
        case 1: return nd::cartesian_product(xc, yv) | vec2 | nd::to_shared();
    }
    throw std::invalid_argument("face_coordinates (invalid axis)");
}

nd::shared_array<numeric::array_t<unit_length, 2>, 2>
cell_coordinates(int block_size, unit_length domain_radius, bsp::tree_index_t<2> block)
{
    static auto memoize = memoizer<int, unit_length, bsp::tree_index_t<2>>(
    [] (auto block_size, auto domain_radius, auto block)
    {
        auto [xv, yv] = vertex_positions(block_size, domain_radius, block);
        auto xc = xv | nd::adjacent_mean(0);
        auto yc = yv | nd::adjacent_mean(0);
        auto vec2 = nd::map(util::apply_to([] (auto x, auto y) { return numeric::array(x, y); }));

        return nd::cartesian_product(xc, yc) | vec2 | nd::to_shared();
    });
    return memoize(block_size, domain_radius, block);
}




//=============================================================================
nd::shared_array<unit_rate, 2>
buffer_rate_field_array(solver_data_t solver_data, bsp::tree_index_t<2> block)
{
    static auto memoize = memoizer<int, unit_length, unit_length, unit_rate, bsp::tree_index_t<2>>(
    [] (auto block_size, auto domain_radius, auto buffer_scale, auto buffer_rate, auto block)
    {
        auto xc = cell_coordinates(block_size, domain_radius, block);
        auto br = xc | nd::map(buffer_rate_field(domain_radius, buffer_scale, buffer_rate));
        return br | nd::to_shared();
    });
    return memoize(solver_data.block_size, solver_data.domain_radius, solver_data.buffer_scale, solver_data.buffer_rate, block);
}

nd::shared_array<iso2d::conserved_density_t, 2>
initial_conserved_array(solver_data_t solver_data, bsp::tree_index_t<2> block)
{
    static auto memoize = memoizer<int, unit_length, unit_length, bsp::tree_index_t<2>>(
    [] (auto block_size, auto domain_radius, auto softening_length, auto block)
    {
        auto xc = cell_coordinates(block_size, domain_radius, block);
        auto uc = xc | nd::map(initial_primitive(softening_length)) | nd::map(iso2d::conserved_density);
        return uc | nd::to_shared();
    });
    return memoize(solver_data.block_size, solver_data.domain_radius, solver_data.softening_length, block);
}

nd::shared_array<numeric::array_t<unit_acceleration, 2>, 2>
centrifugal_term_array(solver_data_t solver_data, bsp::tree_index_t<2> block)
{
    static auto memoize = memoizer<int, unit_length, bsp::tree_index_t<2>>(
    [] (auto block_size, auto domain_radius, auto block)
    {
        auto xc = cell_coordinates(block_size, domain_radius, block);
        auto cen = xc | nd::map(centrifugal_term);
        return cen | nd::to_shared();
    });
    return memoize(solver_data.block_size, solver_data.domain_radius, block);
}




//=============================================================================
solution_t initial(solver_data_t solver_data)
{
    auto uc = bsp::uniform_quadtree(solver_data.depth) |  bsp::maps(std::bind(initial_conserved_array, solver_data, _1));
    return {0, 0.0, uc};
}




//=============================================================================
template<typename ArrayType>
auto extend_x(bsp::shared_tree<pr::computable<ArrayType>, 4> __pc__, bsp::tree_index_t<2> block)
{
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
auto extend_y(bsp::shared_tree<pr::computable<ArrayType>, 4> __pc__, bsp::tree_index_t<2> block)
{
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

auto recover_primitive_array(conserved_array_t uc)
{
    return uc | nd::maps(iso2d::recover_primitive);
}

auto estimate_gradient(primitive_array_t pc, bsp::uint axis, double theta)
{
    return pc | nd::adjacent_zip3(axis) | nd::map(mara::plm_gradient(theta)) | nd::to_shared();
}

auto godunov_fluxes(
    primitive_array_t pc,
    primitive_array_t gc,
    solver_data_t solver_data,
    bsp::tree_index_t<2> block,
    bsp::uint axis)
{
    auto riemann = [axis, rs=solver_data.softening_length] (auto pl, auto pr, auto xf)
    {
        auto cs2 = sound_speed_squared(xf, rs);
        return iso2d::riemann_hlle(pl, pr, cs2, geometric::unit_vector_on(axis + 1));
    };

    auto xf = face_coordinates(solver_data.block_size, solver_data.domain_radius, block, axis);
    auto pl = pc - 0.5 * gc;
    auto pr = pc + 0.5 * gc;

    auto fx = nd::zip(
        pr | nd::select(axis, 0, -1),
        pl | nd::select(axis, 1),
        xf)
    | nd::map(util::apply_to(riemann))
    | nd::to_shared();

    return fx;
}

auto updated_conserved(
    conserved_array_t uc,
    primitive_array_t pc,
    godunov_f_array_t fx,
    godunov_f_array_t fy,
    unit_time time,
    unit_time dt,
    bsp::tree_index_t<2> block,
    solver_data_t solver_data)
{
    auto gravitational_acceleration = [rs=solver_data.softening_length] (auto component)
    {
        return [c=component, rs] (auto p)
        {
            return two_body::gravitational_acceleration(c, p[0], p[1], rs);
        };
    };

    auto xc = cell_coordinates(solver_data.block_size, solver_data.domain_radius, block);
    auto dl = 2.0 * solver_data.domain_radius / double(solver_data.block_size) / double(1 << block.level);

    auto dfx = fx | nd::adjacent_diff(0) | nd::to_shared();
    auto dfy = fy | nd::adjacent_diff(1) | nd::to_shared();

    auto a = unit_length(1.0);
    auto M = unit_mass(1.0);
    auto q = solver_data.mass_ratio;
    auto e = solver_data.eccentricity;
    auto elements = two_body::orbital_elements(a, M, q, e);
    auto inertial = two_body::orbital_state(elements, time);
    auto state    = two_body::rotate(inertial, omega_frame * time);

    auto ag1 = xc | nd::map(gravitational_acceleration(state.first));
    auto ag2 = xc | nd::map(gravitational_acceleration(state.second));
    auto cor = pc | nd::map(iso2d::velocity_vector) | nd::map(coriolis_term);
    auto cen = centrifugal_term_array(solver_data, block);

    auto ss1 = -uc * (xc | nd::map(sink_rate_field(solver_data, two_body::position(state.first))))  | nd::to_shared();
    auto ss2 = -uc * (xc | nd::map(sink_rate_field(solver_data, two_body::position(state.second)))) | nd::to_shared();
    auto sg1 = nd::zip(pc, ag1) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();
    auto sg2 = nd::zip(pc, ag2) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();
    auto sf1 = nd::zip(pc, cen) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();
    auto sf2 = nd::zip(pc, cor) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();

    auto uinit = initial_conserved_array(solver_data, block);
    auto br    = buffer_rate_field_array(solver_data, block);
    auto sb    = -(uc - uinit) * br;
    auto u1    = uc - (dfx + dfy) * dt / dl + (sg1 + sg2 + ss1 + ss2 + sf1 + sf2 + sb) * dt;

    return u1 | nd::to_shared();
}




//=============================================================================
auto encode_step(rational::number_t iteration, unit_time time, bsp::shared_tree<pr::computable<conserved_array_t>, 4> conserved, solver_data_t solver_data)
{
    auto mesh   = indexes(conserved);
    auto _uc_   = conserved;
    auto _pc_   = _uc_ | bsp::maps(pr::map(recover_primitive_array));
    auto _p_ext_x_ = mesh | bsp::maps([_pc_] (auto index) { return extend_x(_pc_, index); });
    auto _p_ext_y_ = mesh | bsp::maps([_pc_] (auto index) { return extend_y(_pc_, index); });
    auto _gradient_x_ = _p_ext_x_ | bsp::maps(pr::map(std::bind(estimate_gradient, _1, 0, 2.0)));
    auto _gradient_y_ = _p_ext_y_ | bsp::maps(pr::map(std::bind(estimate_gradient, _1, 1, 2.0)));
    auto _gradient_ext_x_ = mesh | bsp::maps([_gradient_x_] (auto index) { return extend_x(_gradient_x_, index); });
    auto _gradient_ext_y_ = mesh | bsp::maps([_gradient_y_] (auto index) { return extend_y(_gradient_y_, index); });

    auto _godunov_fluxes_x_ = zip(_p_ext_x_, _gradient_ext_x_, mesh) | bsp::mapvs([solver_data] (auto pc, auto gc, auto block)
    {
        return pr::zip(pc, gc) | pr::mapv(std::bind(godunov_fluxes, _1, _2, solver_data, block, 0));
    });

    auto _godunov_fluxes_y_ = zip(_p_ext_y_, _gradient_ext_y_, mesh) | bsp::mapvs([solver_data] (auto pc, auto gc, auto block)
    {
        return pr::zip(pc, gc) | pr::mapv(std::bind(godunov_fluxes, _1, _2, solver_data, block, 1));
    });

    auto dl = 2.0 * solver_data.domain_radius / double(solver_data.block_size) / double(1 << solver_data.depth);
    auto dt = solver_data.cfl * 0.5 * dl / sqrt(two_body::G * unit_mass(0.5) / solver_data.softening_length);

    auto u1 = zip(_uc_, _pc_, _godunov_fluxes_x_, _godunov_fluxes_y_, mesh)
    | bsp::mapvs([t=time, dt, solver_data] (auto uc, auto pc, auto fx, auto fy, auto block)
    {
        return pr::zip(uc, pc, fx, fy) | pr::mapv(std::bind(updated_conserved, _1, _2, _3, _4, t, dt, block, solver_data));
    });

    return std::tuple(iteration + 1, time + dt, u1);
}




//=============================================================================
solution_t step(solution_t solution, solver_data_t solver_data, int nfold, mara::ThreadPool& pool)
{
    auto [iter, time, uc] = std::tuple(solution.iteration, solution.time, solution.conserved | bsp::maps([] (auto uc) { return pr::just(uc); }));

    for (std::size_t i = 0; i < nfold; ++i)
    {
        std::tie(iter, time, uc) = encode_step(iter, time, uc, solver_data);
    }

    auto master = computable(uc);
    compute(master, pool.scheduler());

    return {
        iter, time, master.value(),
    };
}




//=============================================================================
void side_effects(const mara::config_t& cfg, solution_t solution)
{
    auto cpi = cfg.get_int("cpi");

    if (long(solution.iteration) % cpi == 0)
    {
        auto outdir = cfg.get_string("outdir");
        auto fname = util::format("%s/chkpt.%04ld.h5", outdir.data(), long(solution.iteration) / cpi);

        mara::filesystem::require_dir(outdir);

        auto h5f = h5::File(fname, "w");
        h5::write(h5f, "solution", solution);
        std::printf("write %s\n", fname.data());
    }
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto cfg = config_template().create().update(mara::argv_to_string_map(argc, argv));

    auto solver_data = make_solver_data(cfg);
    auto solution    = initial(solver_data);
    auto orbits      = cfg.get_double("orbits");

    mara::pretty_print(std::cout, "config", cfg);

    auto pool = mara::ThreadPool(cfg.get_int("threads"));
    auto fold = cfg.get_int("fold");

    while (solution.time < orbits * binary_period)
    {
        side_effects(cfg, solution);

        auto [s1, ticks] = control::invoke_timed(step, solution, solver_data, fold, pool);

        solution = s1;

        auto ncells = std::pow(solver_data.block_size * (1 << solver_data.depth), 2);
        auto kzps = 1e6 * fold * ncells / std::chrono::duration_cast<std::chrono::nanoseconds>(ticks).count();

        std::printf("[%06ld] orbit=%lf kzps=%lf\n", long(solution.iteration), solution.time.value / 2 / M_PI, kzps);
        std::fflush(stdout);
    }

    side_effects(cfg, solution);

    return 0;
}
