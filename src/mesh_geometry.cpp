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
#include "mesh_geometry.hpp"




//=============================================================================
mesh::cartesian_1d::geometry_t::geometry_t(dimensional::unit_length domain_size, int block_size)
: domain_size(domain_size)
, block_size(block_size)
{
}

nd::shared_array<mesh::cartesian_1d::coords_t, 1> mesh::cartesian_1d::geometry_t::vert_coordinates(mesh::block_index_t<1> block) const
{
    auto [p0, p1] = block_extent(block);
    auto xv = nd::linspace(p0, p1, block_size + 1);
    return xv | nd::to_shared();
}

nd::shared_array<mesh::cartesian_1d::coords_t, 1> mesh::cartesian_1d::geometry_t::face_coordinates(mesh::block_index_t<1> block, unsigned long axis) const
{
    return invoke_memoized([] (auto self, mesh::block_index_t<1> block, unsigned long axis)
    {
        return std::apply([] (auto... args) { return geometry_t(args...); }, self)
        . vert_coordinates(block)
        | nd::adjacent_mean(axis)
        | nd::to_shared();
    }, as_tuple(), block, axis);
}

nd::shared_array<mesh::cartesian_1d::coords_t, 1> mesh::cartesian_1d::geometry_t::cell_coordinates(mesh::block_index_t<1> block) const
{
    return invoke_memoized([] (auto self, mesh::block_index_t<1> block)
    {
        return std::apply([] (auto... args) { return geometry_t(args...); }, self)
        . vert_coordinates(block)
        | nd::adjacent_mean(0)
        | nd::to_shared();
    }, as_tuple(), block);
}

std::pair<mesh::cartesian_1d::coords_t, mesh::cartesian_1d::coords_t> mesh::cartesian_1d::geometry_t::block_extent(mesh::block_index_t<1> block) const
{
    auto [i0] = numeric::as_tuple(block.coordinates);
    auto [i1] = std::tuple(i0 + 1);
    auto dl = double(1 << block.level);

    auto x0 = domain_size * (-0.5 + 1.0 * i0 / dl);
    auto x1 = domain_size * (-0.5 + 1.0 * i1 / dl);

    return std::pair(x0, x1);
}

mesh::cartesian_1d::coords_t mesh::cartesian_1d::geometry_t::block_centroid(mesh::block_index_t<1> block) const
{
    auto [p0, p1] = block_extent(block);
    return 0.5 * (p0 + p1);
}

dimensional::unit_length mesh::cartesian_1d::geometry_t::cell_spacing(mesh::block_index_t<1> block) const
{
    return domain_size / double(block_size) / double(1 << block.level);
}

std::size_t mesh::cartesian_1d::geometry_t::cells_per_block() const
{
    return block_size;
}

std::tuple<dimensional::unit_length, int> mesh::cartesian_1d::geometry_t::as_tuple() const
{
    return std::tuple(domain_size, block_size);
}




//=============================================================================
mesh::cartesian_2d::geometry_t::geometry_t(dimensional::unit_length domain_size, int block_size)
: domain_size(domain_size)
, block_size(block_size)
{
}

nd::shared_array<mesh::cartesian_2d::coords_t, 2> mesh::cartesian_2d::geometry_t::vert_coordinates(mesh::block_index_t<2> block) const
{
    auto [p0, p1] = block_extent(block);
    auto xv = nd::linspace(p0[0], p1[0], block_size + 1);
    auto yv = nd::linspace(p0[1], p1[1], block_size + 1);
    auto vec2 = nd::mapv([] (auto x, auto y) { return numeric::array(x, y); });

    return nd::cartesian_product(xv, yv) | vec2 | nd::to_shared();
}

nd::shared_array<mesh::cartesian_2d::coords_t, 2> mesh::cartesian_2d::geometry_t::face_coordinates(mesh::block_index_t<2> block, unsigned long axis) const
{
    return invoke_memoized([] (auto self, mesh::block_index_t<2> block, unsigned long axis)
    {
        return std::apply([] (auto... args) { return geometry_t(args...); }, self)
        . vert_coordinates(block)
        | nd::adjacent_mean(1 - axis)
        | nd::to_shared();
    }, as_tuple(), block, axis);
}

nd::shared_array<mesh::cartesian_2d::coords_t, 2> mesh::cartesian_2d::geometry_t::cell_coordinates(mesh::block_index_t<2> block) const
{
    return invoke_memoized([] (auto self, mesh::block_index_t<2> block)
    {
        return std::apply([] (auto... args) { return geometry_t(args...); }, self)
        . vert_coordinates(block)
        | nd::adjacent_mean(0)
        | nd::adjacent_mean(1)
        | nd::to_shared();
    }, as_tuple(), block);
}

std::pair<mesh::cartesian_2d::coords_t, mesh::cartesian_2d::coords_t> mesh::cartesian_2d::geometry_t::block_extent(mesh::block_index_t<2> block) const
{
    auto [i0, j0] = numeric::as_tuple(block.coordinates);
    auto [i1, j1] = std::tuple(i0 + 1, j0 + 1);
    auto dl = double(1 << block.level);

    auto x0 = domain_size * (-0.5 + 1.0 * i0 / dl);
    auto x1 = domain_size * (-0.5 + 1.0 * i1 / dl);
    auto y0 = domain_size * (-0.5 + 1.0 * j0 / dl);
    auto y1 = domain_size * (-0.5 + 1.0 * j1 / dl);

    return std::pair(mesh::cartesian_2d::coords_t{x0, y0}, mesh::cartesian_2d::coords_t{x1, y1});
}

mesh::cartesian_2d::coords_t mesh::cartesian_2d::geometry_t::block_centroid(mesh::block_index_t<2> block) const
{
    auto [p0, p1] = block_extent(block);
    return 0.5 * (p0 + p1);
}

dimensional::unit_length mesh::cartesian_2d::geometry_t::cell_spacing(mesh::block_index_t<2> block) const
{
    return domain_size / double(block_size) / double(1 << block.level);
}

std::size_t mesh::cartesian_2d::geometry_t::cells_per_block() const
{
    return block_size * block_size;
}

std::tuple<dimensional::unit_length, int> mesh::cartesian_2d::geometry_t::as_tuple() const
{
    return std::tuple(domain_size, block_size);
}
