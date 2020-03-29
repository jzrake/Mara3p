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




#include "core_memoize.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"
#include "module_euler.hpp"
#include "scheme_plm_gradient.hpp"




using namespace std::placeholders;
using namespace modules;
using namespace euler2d;




//=============================================================================
template<typename ArrayType>
auto extend_block(bsp::shared_tree<mpr::computable<ArrayType>, 4> tree, bsp::tree_index_t<2> block)
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
            auto ii = i < 2 ? i - 2 + nx : (i >= nx + 2 ? i - 2 - nx : i - 2);
            auto jj = j < 2 ? j - 2 + ny : (j >= ny + 2 ? j - 2 - ny : j - 2);
            return cs[bi][bj](ii, jj);
        }), nd::uivec(nx + 4, ny + 4)) | nd::to_shared();
    });
}




//=============================================================================
mesh_geometry_t::mesh_geometry_t(unit_length domain_size, int block_size)
: domain_size(domain_size)
, block_size(block_size)
{
}

nd::shared_array<coords_t, 2> mesh_geometry_t::vert_coordinates(bsp::tree_index_t<2> block) const
{
    return invoke_memoized([this] (const mesh_geometry_t *self, bsp::tree_index_t<2> block)
    {
        auto [i0, j0] = as_tuple(block.coordinates);
        auto [i1, j1] = std::tuple(i0 + 1, j0 + 1);
        auto dl = double(1 << block.level);

        auto x0 = self->domain_size * (-0.5 + 1.0 * i0 / dl);
        auto x1 = self->domain_size * (-0.5 + 1.0 * i1 / dl);
        auto y0 = self->domain_size * (-0.5 + 1.0 * j0 / dl);
        auto y1 = self->domain_size * (-0.5 + 1.0 * j1 / dl);

        auto xv = nd::linspace(x0, x1, block_size + 1);
        auto yv = nd::linspace(y0, y1, block_size + 1);
        auto vec2 = nd::map(util::apply_to([] (auto x, auto y) { return numeric::array(x, y); }));

        return nd::cartesian_product(xv, yv) | vec2 | nd::to_shared();
    }, this, block);
}

nd::shared_array<coords_t, 2> mesh_geometry_t::face_coordinates(bsp::tree_index_t<2> block, bsp::uint axis) const
{
    return invoke_memoized([] (const mesh_geometry_t *self, bsp::tree_index_t<2> block, bsp::uint axis)
    {
        return self->vert_coordinates(block)
        | nd::adjacent_mean(axis)
        | nd::to_shared();
    }, this, block, axis);
}

nd::shared_array<coords_t, 2> mesh_geometry_t::cell_coordinates(bsp::tree_index_t<2> block) const
{
    return invoke_memoized([] (const mesh_geometry_t *self, bsp::tree_index_t<2> block)
    {
        return self->vert_coordinates(block)
        | nd::adjacent_mean(0)
        | nd::adjacent_mean(1)
        | nd::to_shared();
    }, this, block);
}

unit_length mesh_geometry_t::cell_size(bsp::tree_index_t<2> block) const
{
    return domain_size / double(block_size) / double(1 << block.level);
}




//=============================================================================
conserved_array_t euler2d::initial_conserved_array(
    mesh_geometry_t mesh_geometry,
    bsp::tree_index_t<2> block,
    primitive_mapping_t initial,
    double gamma_law_index)
{
    return mesh_geometry.cell_coordinates(block)
    | nd::map(initial)
    | nd::map(std::bind(euler::conserved_density, _1, gamma_law_index))
    | nd::to_shared();
}




//=============================================================================
conserved_tree_t euler2d::initial_conserved_tree(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry,
    primitive_mapping_t initial,
    double gamma_law_index)
{
    auto u0 = std::bind(initial_conserved_array, mesh_geometry, _1, initial, gamma_law_index);
    auto uc = [u0] (auto block) { return mpr::from(std::bind(u0, block)).name("U"); };
    return mesh_topology | bsp::maps(uc);
}




//=============================================================================
primitive_array_t euler2d::recover_primitive_array(conserved_array_t uc, double gamma_law_index)
{
    return uc
    | nd::map(std::bind(euler::recover_primitive, _1, gamma_law_index))
    | nd::to_shared();
}




//=============================================================================
primitive_array_t euler2d::estimate_gradient(primitive_array_t pc, bsp::uint axis, double plm_theta)
{
    return pc
    | nd::adjacent_zip3(axis)
    | nd::map(mara::plm_gradient(plm_theta))
    | nd::to_shared();
}




//=============================================================================
conserved_array_t euler2d::updated_conserved(
    conserved_array_t uc,
    primitive_array_t pe,
    unit_time time,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    bsp::tree_index_t<2> block,
    double plm_theta,
    double gamma_law_index)
{
    auto dl = mesh_geometry.cell_size(block);
    auto gx = estimate_gradient(pe, 0, plm_theta);
    auto gy = estimate_gradient(pe, 1, plm_theta);
    auto pc = pe | nd::select(0, 2, -2) | nd::select(1, 2, -2) | nd::to_shared();

    auto pe_x = pe | nd::select(1, 2, -2) | nd::select(0, 1, -1);
    auto gx_x = gx | nd::select(1, 2, -2);
    auto pe_y = pe | nd::select(0, 2, -2) | nd::select(1, 1, -1);
    auto gy_y = gy | nd::select(0, 2, -2);

    auto riemann_x = std::bind(euler::riemann_hlle, _1, _2, geometric::unit_vector_on(1), gamma_law_index);
    auto riemann_y = std::bind(euler::riemann_hlle, _1, _2, geometric::unit_vector_on(2), gamma_law_index);

    auto pl_x = (pe_x - 0.5 * gx_x) | nd::select(0, 1);
    auto pr_x = (pe_x + 0.5 * gx_x) | nd::select(0, 0, -1);
    auto pl_y = (pe_y - 0.5 * gy_y) | nd::select(1, 1);
    auto pr_y = (pe_y + 0.5 * gy_y) | nd::select(1, 0, -1);

    auto fx = nd::zip(pr_x, pl_x) | nd::map(util::apply_to(riemann_x));
    auto fy = nd::zip(pr_y, pl_y) | nd::map(util::apply_to(riemann_y));

    auto dfx = fx | nd::adjacent_diff(0);
    auto dfy = fy | nd::adjacent_diff(1);

    return (uc - (dfx + dfy) * dt / dl) | nd::to_shared();
}




//=============================================================================
solution_t euler2d::updated_solution(
    solution_t solution,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    double plm_theta,
    double gamma_law_index)
{
    auto F = std::bind(updated_conserved, _1, _2, solution.time, dt, mesh_geometry, _3, plm_theta, gamma_law_index);
    auto U = [F] (auto b) { return std::bind(F, _1, _2, b); };

    auto mesh = indexes(solution.conserved);
    auto uc   = solution.conserved;
    auto pc   = uc   | bsp::maps(mpr::map([g=gamma_law_index] (auto u) { return recover_primitive_array(u, g); }, "P"));
    auto pe   = mesh | bsp::maps([pc] (auto b) { return extend_block(pc, b); });
    auto u1   = mesh | bsp::maps([uc, pe, U] (auto b) { return zip(value_at(uc, b), value_at(pe, b)) | mpr::mapv(U(b)); });

    return solution_t{solution.iteration + 1, solution.time + dt, u1};
}




//=============================================================================
unit_length euler2d::smallest_cell_size(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry)
{
    return reduce(mesh_topology, [mesh_geometry] (auto seed, auto block)
    {
        return std::min(seed, mesh_geometry.cell_size(block));
    }, unit_length(std::numeric_limits<double>::max()));
}
