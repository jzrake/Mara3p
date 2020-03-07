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




#include <mutex>
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"
#include "physics_two_body.hpp"
#include "scheme_plm_gradient.hpp"
#include "minidisk.hpp"




//=============================================================================
using namespace std::placeholders;
using namespace minidisk;




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
static auto initial_primitive(unit_length softening_length, unit_rate omega_frame)
{
    return [rs=softening_length, omega_frame] (numeric::array_t<unit_length, 2> p)
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

static auto vertex_positions(int block_size, unit_length domain_radius, bsp::tree_index_t<2> block)
{
    auto [i0, j0] = as_tuple(block.coordinates);
    auto [i1, j1] = std::tuple(i0 + 1, j0 + 1);
    auto dl = double(1 << block.level);

    auto x0 = domain_radius * (-1 + 2 * i0 / dl);
    auto x1 = domain_radius * (-1 + 2 * i1 / dl);
    auto y0 = domain_radius * (-1 + 2 * j0 / dl);
    auto y1 = domain_radius * (-1 + 2 * j1 / dl);

    auto xv = nd::linspace(x0, x1, block_size + 1);
    auto yv = nd::linspace(y0, y1, block_size + 1);

    return std::tuple(xv, yv);
}




//=============================================================================
unit_specific_energy minidisk::sound_speed_squared(numeric::array_t<unit_length, 2> p, unit_length softening_length, unit_scalar mach_number)
{
    auto [x, y] = as_tuple(p);
    auto ph = two_body::potential(two_body::point_mass_t(), x, y, softening_length);
    return -ph / mach_number / mach_number;
}

nd::shared_array<numeric::array_t<unit_length, 2>, 2>
minidisk::face_coordinates(int block_size, unit_length domain_radius, bsp::tree_index_t<2> block, bsp::uint axis)
{
    static auto memoize = memoizer<int, unit_length, bsp::tree_index_t<2>, bsp::uint>([] (auto block_size, auto domain_radius, auto block, auto axis)
    {
        auto [xv, yv] = vertex_positions(block_size, domain_radius, block);
        auto xc = xv | nd::adjacent_mean(0);
        auto yc = yv | nd::adjacent_mean(0);
        auto vec2 = nd::map(util::apply_to([] (auto x, auto y) { return numeric::array(x, y); }));

        switch (axis)
        {
            case 0: return nd::cartesian_product(xv, yc) | vec2 | nd::to_shared();
            case 1: return nd::cartesian_product(xc, yv) | vec2 | nd::to_shared();
        }
        throw std::invalid_argument("face_coordinates (invalid axis)");
    });
    return memoize(block_size, domain_radius, block, axis);
}

nd::shared_array<numeric::array_t<unit_length, 2>, 2>
minidisk::cell_coordinates(int block_size, unit_length domain_radius, bsp::tree_index_t<2> block)
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
minidisk::buffer_rate_field_array(solver_data_t solver_data, bsp::tree_index_t<2> block)
{
    static auto memoize = memoizer<int, unit_length, unit_length, unit_rate, bsp::tree_index_t<2>>(
    [] (auto block_size, auto domain_radius, auto buffer_scale, auto buffer_rate, auto block)
    {
        auto xc = cell_coordinates(block_size, domain_radius, block);
        auto br = xc | nd::map(buffer_rate_field(domain_radius, buffer_scale, buffer_rate));
        return br | nd::to_shared();
    });
    return memoize(
        solver_data.block_size,
        solver_data.domain_radius,
        solver_data.buffer_scale,
        solver_data.buffer_rate, block);
}

nd::shared_array<iso2d::conserved_density_t, 2>
minidisk::initial_conserved_array(solver_data_t solver_data, bsp::tree_index_t<2> block)
{
    static auto memoize = memoizer<int, unit_length, unit_length, unit_rate, bsp::tree_index_t<2>>(
    [] (auto block_size, auto domain_radius, auto softening_length, auto omega_frame, auto block)
    {
        auto xc = cell_coordinates(block_size, domain_radius, block);
        auto uc = xc | nd::map(initial_primitive(softening_length, omega_frame)) | nd::map(iso2d::conserved_density);
        return uc | nd::to_shared();
    });
    return memoize(
        solver_data.block_size,
        solver_data.domain_radius,
        solver_data.softening_length,
        solver_data.omega_frame,
        block);
}

nd::shared_array<numeric::array_t<unit_acceleration, 2>, 2>
minidisk::centrifugal_term_array(solver_data_t solver_data, bsp::tree_index_t<2> block)
{
    static auto memoize = memoizer<int, unit_length, unit_rate, bsp::tree_index_t<2>>(
    [] (auto block_size, auto domain_radius, auto omega_frame, auto block)
    {
        auto xc = cell_coordinates(block_size, domain_radius, block);
        auto cen = xc | nd::map(std::bind(centrifugal_term, _1, omega_frame));
        return cen | nd::to_shared();
    });
    return memoize(
        solver_data.block_size,
        solver_data.domain_radius,
        solver_data.omega_frame,
        block);
}




//=============================================================================
primitive_array_t minidisk::recover_primitive_array(conserved_array_t uc)
{
    return uc | nd::maps(iso2d::recover_primitive);
}

primitive_array_t minidisk::estimate_gradient(primitive_array_t pc, bsp::uint axis, double theta)
{
    return pc | nd::adjacent_zip3(axis) | nd::map(mara::plm_gradient(theta)) | nd::to_shared();
}

godunov_f_array_t minidisk::godunov_fluxes(
    primitive_array_t pc,
    primitive_array_t gc,
    solver_data_t solver_data,
    bsp::tree_index_t<2> block,
    bsp::uint axis)
{
    auto riemann = [axis, rs=solver_data.softening_length, mach=solver_data.mach_number] (auto pl, auto pr, auto xf)
    {
        auto cs2 = sound_speed_squared(xf, rs, mach);
        return iso2d::riemann_hlle(pl, pr, cs2, geometric::unit_vector_on(axis + 1));
    };

    auto xf = minidisk::face_coordinates(solver_data.block_size, solver_data.domain_radius, block, axis);
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

godunov_f_array_t minidisk::godunov_and_viscous_fluxes(
    primitive_array_t pc,
    primitive_array_t gc_long,
    primitive_array_t gc_tran,
    solver_data_t solver_data,
    bsp::tree_index_t<2> block,
    bsp::uint axis)
{
    auto riemann = [axis, rs=solver_data.softening_length, mach=solver_data.mach_number] (auto pl, auto pr, auto xf)
    {
        auto cs2 = sound_speed_squared(xf, rs, mach);
        return iso2d::riemann_hlle(pl, pr, cs2, geometric::unit_vector_on(axis + 1));
    };

    auto xf = face_coordinates(solver_data.block_size, solver_data.domain_radius, block, axis);
    auto pl = (pc - 0.5 * gc_long) | nd::select(axis, 1)     | nd::to_dynamic();
    auto pr = (pc + 0.5 * gc_long) | nd::select(axis, 0, -1) | nd::to_dynamic();
    auto fx = nd::zip(pr, pl, xf) | nd::map(util::apply_to(riemann));

    auto viscous_flux = [=] ()
    {
        auto nu = solver_data.kinematic_viscosity;
        auto sl = map(pl, iso2d::mass_density);
        auto sr = map(pr, iso2d::mass_density);
        auto mu = nu * 0.5 * (sl + sr);
        auto dl = cell_size(block, solver_data);
        auto block_size = solver_data.block_size;

        switch (axis)
        {
            case 0:
            {
                auto dx_vx = gc_long | nd::map(iso2d::velocity_1) | nd::adjacent_mean(0) | nd::divide(dl);
                auto dy_vx = gc_tran | nd::map(iso2d::velocity_1) | nd::adjacent_mean(0) | nd::divide(dl);
                auto dx_vy = gc_long | nd::map(iso2d::velocity_2) | nd::adjacent_mean(0) | nd::divide(dl);
                auto dy_vy = gc_tran | nd::map(iso2d::velocity_2) | nd::adjacent_mean(0) | nd::divide(dl);

                auto tauxs = nd::zeros<dimensional::quantity_t<1,-1,-1>>(block_size + 1, block_size);
                auto tauxx = mu * (dx_vx - dy_vy);
                auto tauxy = mu * (dx_vy + dy_vx);

                return nd::zip(tauxs, -tauxx, -tauxy) | nd::construct<iso2d::flux_vector_t>() | nd::to_dynamic();
            }
            case 1:
            {
                auto dx_vx = gc_tran | nd::map(iso2d::velocity_1) | nd::adjacent_mean(1) | nd::divide(dl);
                auto dy_vx = gc_long | nd::map(iso2d::velocity_1) | nd::adjacent_mean(1) | nd::divide(dl);
                auto dx_vy = gc_tran | nd::map(iso2d::velocity_2) | nd::adjacent_mean(1) | nd::divide(dl);
                auto dy_vy = gc_long | nd::map(iso2d::velocity_2) | nd::adjacent_mean(1) | nd::divide(dl);

                auto tauys = nd::zeros<dimensional::quantity_t<1,-1,-1>>(block_size, block_size + 1);
                auto tauyx = mu * (dx_vy + dy_vx);
                auto tauyy =-mu * (dx_vx - dy_vy);

                return nd::zip(tauys, -tauyx, -tauyy) | nd::construct<iso2d::flux_vector_t>() | nd::to_dynamic();
            }
        }
        throw;
    };

    if (solver_data.kinematic_viscosity == unit_viscosity(0.0))
    {
        return fx | nd::to_shared();
    }
    return (fx + viscous_flux()) | nd::to_shared();
}

conserved_array_t minidisk::updated_conserved(
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
    auto dl = cell_size(block, solver_data);

    auto dfx = fx | nd::adjacent_diff(0);
    auto dfy = fy | nd::adjacent_diff(1);

    auto a = unit_length(1.0);
    auto M = unit_mass(1.0);
    auto q = solver_data.mass_ratio;
    auto e = solver_data.eccentricity;
    auto elements = two_body::orbital_elements(a, M, q, e);
    auto inertial = two_body::orbital_state(elements, time);
    auto state    = two_body::rotate(inertial, solver_data.omega_frame * time);

    auto ag1 = xc | nd::map(gravitational_acceleration(state.first));
    auto ag2 = xc | nd::map(gravitational_acceleration(state.second));
    auto cor = pc | nd::map(iso2d::velocity_vector) | nd::map(std::bind(coriolis_term, _1, solver_data.omega_frame));
    auto cen = centrifugal_term_array(solver_data, block);

    auto ss1 = -uc * (xc | nd::map(sink_rate_field(solver_data, two_body::position(state.first))));
    auto ss2 = -uc * (xc | nd::map(sink_rate_field(solver_data, two_body::position(state.second))));
    auto sg1 = nd::zip(pc, ag1) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));
    auto sg2 = nd::zip(pc, ag2) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));
    auto sf1 = nd::zip(pc, cen) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));
    auto sf2 = nd::zip(pc, cor) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));

    auto uinit = initial_conserved_array(solver_data, block);
    auto br    = buffer_rate_field_array(solver_data, block);
    auto sb    = -(uc - uinit) * br;
    auto u1    = uc - (dfx + dfy) * dt / dl + (sg1 + sg2 + ss1 + ss2 + sf1 + sf2 + sb) * dt;

    return u1 | nd::to_shared();
}
