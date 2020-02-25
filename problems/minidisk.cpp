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




//=============================================================================
using namespace dimensional;
using namespace std::placeholders;

static auto mach_number = 10.0;
static auto omega_frame = dimensional::quantity_t<0, 0, -1>(1.0);




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("block_size",            32)   // number of zones per side
    .item("depth",                  2)   // number of levels in the mesh
    .item("domain_radius",        1.5)   // half-size of the domain
    .item("softening_length",    0.01)   // gravitational softening length
    .item("sink_radius",         0.01)   // radius of mass sinks
    .item("sink_rate",            1e3)   // rate of mass and momentum removal at sinks
    .item("buffer_rate",          1e2)   // maximum rate of buffer driving
    .item("buffer_scale",         0.2)   // buffer onset distance
    .item("tfinal",             100.0)   // time to stop the simulation
    .item("threads",                1)   // number of threads to execute on
    .item("restart", std::string(""))    // a checkpoint file to restart from
    .item("outdir",  std::string("."));  // the directory where output files are written
}

struct solver_data_t
{
    int                 block_size;
    int                 depth;
    unit_length         domain_radius;
    unit_length         softening_length;
    unit_length         sink_radius;
    unit_rate           sink_rate;
    unit_rate           buffer_rate;
    unit_length         buffer_scale;
};

auto make_solver_data(const mara::config_t& cfg)
{
    return solver_data_t{
        cfg.get_int("block_size"),
        cfg.get_int("depth"),
        cfg.get_double("domain_radius"),
        cfg.get_double("softening_length"),
        cfg.get_double("sink_radius"),
        cfg.get_double("sink_rate"),
        cfg.get_double("buffer_rate"),
        cfg.get_double("buffer_scale"),
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
auto initial_primitive(solver_data_t solver_data)
{
    return [rs=solver_data.softening_length] (numeric::array_t<unit_length, 2> p)
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

auto buffer_rate_field(solver_data_t solver_data)
{
    return [solver_data] (numeric::array_t<unit_length, 2> p)
    {
        auto r = sqrt(sum(p * p));
        auto y = (r - solver_data.domain_radius) / solver_data.buffer_scale;
        return 0.5 * solver_data.buffer_rate * (1.0 + std::tanh(y));
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

auto cell_coordinates(bsp::uint n, unit_length dr, bsp::tree_index_t<2> block)
{
    auto [xv, yv] = vertex_positions(n, dr, block);
    auto xc = xv | nd::adjacent_mean(0);
    auto yc = yv | nd::adjacent_mean(0);
    auto vec2 = nd::map(util::apply_to([] (auto x, auto y) { return numeric::array(x, y); }));

    return nd::cartesian_product(xc, yc) | vec2;
}




//=============================================================================
solution_t initial(solver_data_t solver_data)
{
    auto iu = [solver_data] (auto block)
    {
        auto xc = cell_coordinates(solver_data.block_size, solver_data.domain_radius, block);
        auto uc = xc | nd::map(initial_primitive(solver_data)) | nd::map(iso2d::conserved_density) | nd::to_shared();
        return uc;
    };

    auto _tp_ = bsp::uniform_quadtree(solver_data.depth);
    auto _uc_ = _tp_ | bsp::map(iu) | bsp::to_shared();

    return {0, 0.0, _uc_};
}




//=============================================================================
template<typename ArrayType>
auto extend(bsp::shared_tree<pr::computable<ArrayType>, 4> __pc__, bsp::tree_index_t<2> block, bsp::uint axis)
{
    auto _pl_ = value_at(__pc__, prev_on(block, axis));
    auto _pc_ = value_at(__pc__, block);
    auto _pr_ = value_at(__pc__, next_on(block, axis));

    return pr::zip(_pl_, _pc_, _pr_) | pr::mapv([axis] (auto pl, auto pc, auto pr)
    {
        return select(pl, axis, -1) | nd::concat(pc, axis) | nd::concat(select(pr, axis, 0, 1), axis) | nd::to_shared();
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

    auto dfx = fx | nd::adjacent_diff(0);
    auto dfy = fy | nd::adjacent_diff(1);

    auto elements = two_body::orbital_elements();
    auto inertial = two_body::orbital_state(elements, time);
    auto state    = two_body::rotate(inertial, omega_frame * time);

    auto ag1 = xc | nd::map(gravitational_acceleration(state.first));
    auto ag2 = xc | nd::map(gravitational_acceleration(state.second));
    auto cen = xc | nd::map(centrifugal_term);
    auto cor = pc | nd::map(iso2d::velocity_vector) | nd::map(coriolis_term);

    auto ss1 = -uc * (xc | nd::map(sink_rate_field(solver_data, two_body::position(state.first))))  | nd::to_shared();
    auto ss2 = -uc * (xc | nd::map(sink_rate_field(solver_data, two_body::position(state.second)))) | nd::to_shared();
    auto sg1 = nd::zip(pc, ag1) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();
    auto sg2 = nd::zip(pc, ag2) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();
    auto sf1 = nd::zip(pc, cen) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();
    auto sf2 = nd::zip(pc, cor) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms))      | nd::to_shared();

    auto uinit = xc | nd::map(initial_primitive(solver_data)) | nd::map(iso2d::conserved_density);
    auto br    = xc | nd::map(buffer_rate_field(solver_data));
    auto sb    = -(uc - uinit) * br;
    auto u1    = uc - (dfx + dfy) * dt / dl + (sg1 + sg2 + ss1 + ss2 + sf1 + sf2 + sb) * dt;

    return u1 | nd::to_shared();
}




//=============================================================================
void side_effects(solution_t solution)
{
    if (long(solution.iteration) % 50 == 0)
    {
        auto fname = util::format("chkpt.%04ld.h5", long(solution.iteration) / 50);
        auto h5f = h5::File(fname, "w");
        h5::write(h5f, "solution", solution);
        std::printf("write %s\n", fname.data());
    }
}




//=============================================================================
solution_t step(solution_t solution, solver_data_t solver_data, mara::ThreadPool& pool)
{
    auto mesh   = indexes(solution.conserved);
    auto _uc_   = solution.conserved | bsp::maps([] (auto uc) { return pr::just(uc); });
    auto _pc_   = _uc_ | bsp::maps(pr::map(recover_primitive_array));

    auto _p_ext_x_ = mesh | bsp::maps([_pc_] (auto index) { return extend(_pc_, index, 0); });
    auto _p_ext_y_ = mesh | bsp::maps([_pc_] (auto index) { return extend(_pc_, index, 1); });
    auto _gradient_x_ = _p_ext_x_ | bsp::maps(pr::map(std::bind(estimate_gradient, _1, 0, 2.0)));
    auto _gradient_y_ = _p_ext_y_ | bsp::maps(pr::map(std::bind(estimate_gradient, _1, 1, 2.0)));
    auto _gradient_ext_x_ = mesh | bsp::maps([_gradient_x_] (auto index) { return extend(_gradient_x_, index, 0); });
    auto _gradient_ext_y_ = mesh | bsp::maps([_gradient_y_] (auto index) { return extend(_gradient_y_, index, 1); });

    auto _godunov_fluxes_x_ = zip(_p_ext_x_, _gradient_ext_x_, mesh) | bsp::mapvs([solver_data] (auto pc, auto gc, auto i)
    {
        return pr::zip(pc, gc) | pr::mapv(std::bind(godunov_fluxes, _1, _2, solver_data, i, 0));
    });

    auto _godunov_fluxes_y_ = zip(_p_ext_y_, _gradient_ext_y_, mesh) | bsp::mapvs([solver_data] (auto pc, auto gc, auto i)
    {
        return pr::zip(pc, gc) | pr::mapv(std::bind(godunov_fluxes, _1, _2, solver_data, i, 1));
    });

    auto dl = 2.0 * solver_data.domain_radius / double(solver_data.block_size) / double(1 << solver_data.depth);
    auto dt = 0.25 * dl / unit_velocity(10.0);

    auto u1 = computable(zip(_uc_, _pc_, _godunov_fluxes_x_, _godunov_fluxes_y_, mesh)
    | bsp::mapvs([t=solution.time, dt, solver_data] (auto uc, auto pc, auto fx, auto fy, auto block)
    {
        return pr::zip(uc, pc, fx, fy) | pr::mapv(std::bind(updated_conserved, _1, _2, _3, _4, t, dt, block, solver_data));
    }));

    compute(u1, pool.scheduler());

    return {
        solution.iteration + 1,
        solution.time + dt,
        solution.conserved = u1.value(),
    };
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto cfg = config_template().create().update(mara::argv_to_string_map(argc, argv));

    auto solver_data = make_solver_data(cfg);
    auto solution    = initial(solver_data);
    auto tfinal      = unit_time(cfg.get_double("tfinal"));

    mara::pretty_print(std::cout, "config", cfg);

    auto pool = mara::ThreadPool(cfg.get_int("threads"));

    while (solution.time < tfinal)
    {
        side_effects(solution);

        auto [s1, ticks] = control::invoke_timed(step, solution, solver_data, pool);

        solution = s1;

        auto ncells = std::pow(solver_data.block_size * (1 << solver_data.depth), 2);
        auto kzps = 1e6 * ncells / std::chrono::duration_cast<std::chrono::nanoseconds>(ticks).count();

        std::printf("[%06ld] orbit=%lf kzps=%lf\n", long(solution.iteration), solution.time.value / 2 / M_PI, kzps);
    }

    return 0;
}
