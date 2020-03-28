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
#include "scheme_euler.hpp"
#include "scheme_plm_gradient.hpp"




using namespace std::placeholders;
using namespace schemes;
using namespace euler2d;




//=============================================================================
mesh_coordinates_t::mesh_coordinates_t(unit_length domain_size, int block_size)
: domain_size(domain_size)
, block_size(block_size)
{
}

nd::shared_array<coords_t, 2> mesh_coordinates_t::vert_coordinates(bsp::tree_index_t<2> block) const
{
    return invoke_memoized([this] (const mesh_coordinates_t *self, bsp::tree_index_t<2> block)
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

nd::shared_array<coords_t, 2> mesh_coordinates_t::face_coordinates(bsp::tree_index_t<2> block, bsp::uint axis) const
{
    return invoke_memoized([] (const mesh_coordinates_t *self, bsp::tree_index_t<2> block, bsp::uint axis)
    {
        return self->vert_coordinates(block)
        | nd::adjacent_mean(axis)
        | nd::to_shared();
    }, this, block, axis);
}

nd::shared_array<coords_t, 2> mesh_coordinates_t::cell_coordinates(bsp::tree_index_t<2> block) const
{
    return invoke_memoized([] (const mesh_coordinates_t *self, bsp::tree_index_t<2> block)
    {
        return self->vert_coordinates(block)
        | nd::adjacent_mean(0)
        | nd::adjacent_mean(1)
        | nd::to_shared();
    }, this, block);
}

unit_length mesh_coordinates_t::cell_size(bsp::tree_index_t<2> block) const
{
    return domain_size / double(block_size) / double(1 << block.level);
}




//=============================================================================
conserved_array_t scheme_t::initial_conserved_array(bsp::tree_index_t<2> block, primitive_mapping_t initial) const
{
    return mesh.cell_coordinates(block)
    | nd::map(initial)
    | nd::map(std::bind(euler::conserved_density, _1, gamma_law_index))
    | nd::to_shared();
}

primitive_array_t scheme_t::recover_primitive_array(conserved_array_t uc) const
{
    return uc
    | nd::map(std::bind(euler::recover_primitive, _1, gamma_law_index))
    | nd::to_shared();
}

primitive_array_t scheme_t::estimate_gradient(primitive_array_t pc, bsp::uint axis) const
{
    return pc
    | nd::adjacent_zip3(axis)
    | nd::map(mara::plm_gradient(plm_theta))
    | nd::to_shared();
}

conserved_array_t scheme_t::updated_conserved(
    conserved_array_t uc,
    primitive_array_t pe,
    unit_time time,
    unit_time dt,
    bsp::tree_index_t<2> block) const
{
    auto dl = mesh.cell_size(block);
    auto gx = estimate_gradient(pe, 0);
    auto gy = estimate_gradient(pe, 1);
    auto pc = pe | nd::select(0, 2, -2) | nd::select(1, 2, -2) | nd::to_shared();

    auto pe_x = pe | nd::select(1, 2, -2) | nd::select(0, 1, -1);
    auto gx_x = gx | nd::select(1, 2, -2);
    auto gy_x = gy | nd::select(1, 1, -1) | nd::select(0, 1, -1);
    auto pe_y = pe | nd::select(0, 2, -2) | nd::select(1, 1, -1);
    auto gy_y = gy | nd::select(0, 2, -2);
    auto gx_y = gx | nd::select(0, 1, -1) | nd::select(1, 1, -1);

    auto riemann_x = std::bind(euler::riemann_hlle, _1, _2, geometric::unit_vector_on(0), gamma_law_index);
    auto riemann_y = std::bind(euler::riemann_hlle, _1, _2, geometric::unit_vector_on(1), gamma_law_index);

    auto pl_x = (pc - 0.5 * gx_x) | nd::select(0, 1);
    auto pr_x = (pc + 0.5 * gx_x) | nd::select(0, 0, -1);
    auto pl_y = (pc - 0.5 * gy_y) | nd::select(1, 1);
    auto pr_y = (pc + 0.5 * gy_y) | nd::select(1, 0, -1);

    auto fx = nd::zip(pr_x, pl_x) | nd::map(util::apply_to(riemann_x));
    auto fy = nd::zip(pr_y, pl_y) | nd::map(util::apply_to(riemann_y));

    auto dfx = fx | nd::adjacent_diff(0);
    auto dfy = fy | nd::adjacent_diff(1);

    auto u1 = uc - (dfx + dfy) * dt / dl;

    return u1 | nd::to_shared();
}
