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




#include "mesh_block_index.hpp"
#include "core_numeric_array.hpp"
#include "core_dimensional.hpp"
#include "core_ndarray.hpp"




//=============================================================================
namespace mesh
{




//=============================================================================
namespace cartesian_1d
{


using coords_t = dimensional::unit_length;


//=============================================================================
struct geometry_t
{
    geometry_t() {}
    geometry_t(dimensional::unit_length domain_size, int block_size);

    nd::shared_array<coords_t, 1> vert_coordinates(mesh::block_index_t<1> block) const;
    nd::shared_array<coords_t, 1> face_coordinates(mesh::block_index_t<1> block, unsigned long axis) const;
    nd::shared_array<coords_t, 1> cell_coordinates(mesh::block_index_t<1> block) const;
    std::pair<coords_t, coords_t> block_extent(mesh::block_index_t<1> block) const;
    coords_t block_centroid(mesh::block_index_t<1> block) const;
    dimensional::unit_length cell_spacing(mesh::block_index_t<1> block) const;
    std::size_t cells_per_block() const;
    std::tuple<dimensional::unit_length, int> as_tuple() const;

    dimensional::unit_length domain_size = 1.0;
    int block_size = 64;
};

} // cartesian_1d




//=============================================================================
namespace cartesian_2d
{


using coords_t = numeric::array_t<dimensional::unit_length, 2>;


//=============================================================================
struct geometry_t
{
    geometry_t() {}
    geometry_t(dimensional::unit_length domain_size, int block_size);

    nd::shared_array<coords_t, 2> vert_coordinates(mesh::block_index_t<2> block) const;
    nd::shared_array<coords_t, 2> face_coordinates(mesh::block_index_t<2> block, unsigned long axis) const;
    nd::shared_array<coords_t, 2> cell_coordinates(mesh::block_index_t<2> block) const;
    std::pair<coords_t, coords_t> block_extent(mesh::block_index_t<2> block) const;
    coords_t block_centroid(mesh::block_index_t<2> block) const;
    dimensional::unit_length cell_spacing(mesh::block_index_t<2> block) const;
    std::size_t cells_per_block() const;
    std::tuple<dimensional::unit_length, int> as_tuple() const;

    dimensional::unit_length domain_size = 1.0;
    int block_size = 64;
};

} // cartesian_2d

} // namespace mesh2d
