/**
 ==============================================================================
 Copyright 2019, Jonathan Zrake

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




#pragma once
#include <map>
#include <string>
#include <vector>
#include "core_geometric.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_numeric_array.hpp"
#include "core_sequence.hpp"
#include "core_util.hpp"




//=============================================================================
namespace mesh {




/**
 * @brief      This enum is used to represent an axis in a 3d array
 */
enum class axis_3d {
    i, j, k
};




//=============================================================================
inline auto to_numeric_array(nd::uivec_t<3> i)
{
    return numeric::array(i[0], i[1], i[2]);
};

inline auto to_uivec(numeric::array_t<unsigned long, 3> a)
{
    return nd::uivec(a[0], a[1], a[2]);
};

inline auto kronecker_delta(axis_3d axis)
{
    switch (axis)
    {
        case axis_3d::i: return numeric::array(1, 0, 0);
        case axis_3d::j: return numeric::array(0, 1, 0);
        case axis_3d::k: return numeric::array(0, 0, 1);
        default: return numeric::array(0, 0, 0);
    }
}

inline nd::uivec_t<3> block_extent(nd::uint depth)
{
    return nd::uivec(1 << depth, 1 << depth, 1 << depth);
}




/**
 * @brief      Return a 3D lattice of Euclidean vectors of a floating point data
 *             type T (e.g. float, double, dimensional::unit_length), evenly
 *             spaced between 0 and 1 on each axis. These points represent the
 *             location of mesh vertices, so the endpoints are included.
 *
 * @param[in]  ni    The number of lattice points in the x-direction
 * @param[in]  nj    The number of lattice points in the y-direction
 * @param[in]  nk    The number of lattice points in the z-direction
 *
 * @tparam     T     The data type of the Euclidean vectors
 *
 * @return     A 3D array of 3D vectors
 */
template<typename T=double>
auto unit_lattice(nd::uint ni, nd::uint nj, nd::uint nk)
{
    return cartesian_product(
        nd::linspace(0.0, 1.0, ni),
        nd::linspace(0.0, 1.0, nj),
        nd::linspace(0.0, 1.0, nk))
    | nd::map(util::apply_to(geometric::euclidean_vector<T>));
}




/**
 * @brief      Return a 3D array of vertex location vectors at the given index
 *             of a multi-level mesh.
 *
 * @param[in]  level        The block level
 * @param[in]  coordinates  The block coordinates
 * @param[in]  block_size   The number of cells per block, per side
 *
 * @tparam     T            The data type used for the euclidean vectors
 *
 * @return     A 3D array of 3D vectors
 */
template<typename T=double>
auto construct_vertices(unsigned long level, numeric::array_t<unsigned long, 3> coordinates, nd::uint block_size)
{
    auto N = block_size;
    auto dx = T(1.0 / (1 << level));
    auto x0 = dx * geometric::to_euclidean_vector(numeric::construct<double>(coordinates));
    auto xv = dx * unit_lattice(N + 1, N + 1, N + 1) + x0;
    return xv;
}




/**
 * @brief      Return a 3-tuple of edge position arrays, given a 3D array of 3D
 *             vertices. The edge positions are at the geometric midpoint
 *             between two vertices. The array of i-oriented edge positions is
 *             one element smaller on axis i than the vertex array.
 *
 * @param[in]  vertices   The array of vertices
 *
 * @tparam     P          The provider type of the vertex array
 * @tparam     T          The value type of the Euclidean vectors (inferred)
 * @tparam     <unnamed>  Enable if vertex array is of Euclidean vectors
 *
 * @return     A 3-tuple of edge position arrays
 */
template<typename P,
    typename T = typename nd::array_t<P, 3>::value_type,
    typename = std::enable_if_t<geometric::is_euclidean_vector<T>::value>>
auto edge_positions(nd::array_t<P, 3> vertices)
{
    auto a1 = nd::adjacent_mean(0);
    auto a2 = nd::adjacent_mean(1);
    auto a3 = nd::adjacent_mean(2);
    return std::tuple(a1(vertices), a2(vertices), a3(vertices));
}




/**
 * @brief      Return a 3-tuple of face position arrays, given a 3D array of 3D
 *             vertices. The face positions are at the geometric centroid of
 *             four vertices. The array of i-oriented face positions has the
 *             same length on axis i as the vertex array, but is one element
 *             smaller on the other two axes.
 *
 * @param[in]  vertices   The array of vertices
 *
 * @tparam     P          The provider type of the vertex array
 * @tparam     T          The value type of the Euclidean vectors (inferred)
 * @tparam     <unnamed>  Enable if vertex array is of Euclidean vectors
 *
 * @return     A 3-tuple of edge position arrays
 */
template<typename P,
    typename T = typename nd::array_t<P, 3>::value_type,
    typename = std::enable_if_t<geometric::is_euclidean_vector<T>::value>>
auto face_positions(nd::array_t<P, 3> vertices)
{
    auto a1 = nd::adjacent_mean(0);
    auto a2 = nd::adjacent_mean(1);
    auto a3 = nd::adjacent_mean(2);
    return std::tuple(vertices | a2 | a3, vertices | a3 | a1, vertices | a1 | a2);
}




/**
 * @brief      { function_description }
 *
 * @param[in]  vertices   The vertices
 *
 * @tparam     P          { description }
 * @tparam     T          { description }
 * @tparam     <unnamed>  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename P,
    typename T = typename nd::array_t<P, 3>::value_type,
    typename = std::enable_if_t<geometric::is_euclidean_vector<T>::value>>
auto cell_positions(nd::array_t<P, 3> vertices)
{
    auto a1 = nd::adjacent_mean(0);
    auto a2 = nd::adjacent_mean(1);
    auto a3 = nd::adjacent_mean(2);
    return vertices | a1 | a2 | a3;
}




/**
 * @brief      { function_description }
 *
 * @param[in]  b1    The b 1
 * @param[in]  b2    The b 2
 * @param[in]  b3    The b 3
 *
 * @tparam     P     { description }
 * @tparam     T     { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename P, typename T = typename nd::array_t<P, 3>::value_type>
auto face_to_cell(nd::array_t<P, 3> b1, nd::array_t<P, 3> b2, nd::array_t<P, 3> b3)
{
    auto a1 = nd::adjacent_mean(0);
    auto a2 = nd::adjacent_mean(1);
    auto a3 = nd::adjacent_mean(2);
    return nd::zip(b1 | a1, b2 | a2, b3 | a3) | nd::map(util::apply_to(geometric::euclidean_vector<T>));
}




/**
 * @brief      Evaluate the area-integrated curl of a vector field, discretized
 *             on a 3D, logically cartesian staggered mesh. The argument arrays
 *             are the line-integrals E.dl of some vector field E, on the mesh
 *             edges. This function evaluates the integral of curl(E).dA over
 *             the mesh faces, and returns three arrays, for the x, y, and
 *             z-oriented faces.
 *
 * @param[in]  e1    E.dl on the x-oriented edges
 * @param[in]  e2    E.dl on the y-oriented edges
 * @param[in]  e3    E.dl on the z-oriented edges
 *
 * @tparam     T     The provider type of the arrays
 *
 * @return     A tuple of three arrays, containing the surface-integral of
 *             curl(E).dA, computed on the faces of the corresponding logical
 *             axes.
 */
template<typename P>
auto solenoidal_difference(nd::array_t<P, 3> e1, nd::array_t<P, 3> e2, nd::array_t<P, 3> e3)
{
    auto d1 = nd::adjacent_diff(0);
    auto d2 = nd::adjacent_diff(1);
    auto d3 = nd::adjacent_diff(2);
    return std::tuple(d2(e3) - d3(e2), d3(e1) - d1(e3), d1(e2) - d2(e1));
}




/**
 * @brief      Evaluate the volume-integrated divergence of a vector field,
 *             discretized a 3D, logically cartesian mesh. The argument arrays
 *             are the area-integrals F.dA of some vector field F, on the mesh
 *             faces. This function evaluates the integral of div(F) dV over the
 *             mesh volumes.
 *
 * @param[in]  f1    The F.dA on the x-oriented faces
 * @param[in]  f2    The F.dA on the y-oriented faces
 * @param[in]  f3    The F.dA on the z-oriented faces
 *
 * @tparam     P1    The provider type of the array 1
 * @tparam     P2    The provider type of the array 2
 * @tparam     P3    The provider type of the array 3
 *
 * @return     A single 3D array
 */
template<typename P1, typename P2, typename P3>
auto divergence_difference(nd::array_t<P1, 3> f1, nd::array_t<P2, 3> f2, nd::array_t<P3, 3> f3)
{
    auto d1 = nd::adjacent_diff(0);
    auto d2 = nd::adjacent_diff(1);
    auto d3 = nd::adjacent_diff(2);
    return d1(f1) + d2(f2) + d3(f3);
}




//=============================================================================
template<typename P, typename = std::enable_if_t<nd::is_array<std::invoke_result_t<P, nd::uivec_t<3>>>::value>>
auto tile_blocks(nd::array_t<P, 3> blocks)
{
    for (auto block : blocks)
        if (shape(block) != shape(front(blocks)))
            throw std::invalid_argument("mesh::tile_blocks (supplied blocks have non-uniform shapes)");

    auto bi = shape(blocks, 0);
    auto bj = shape(blocks, 1);
    auto bk = shape(blocks, 2);
    auto ni = shape(front(blocks), 0);
    auto nj = shape(front(blocks), 1);
    auto nk = shape(front(blocks), 2);

    return nd::make_array(nd::indexing([ni, nj, nk, blocks] (auto i, auto j, auto k)
    {
        auto b = nd::uivec(i / ni, j / nj, k / nk);
        auto c = nd::uivec(i % ni, j % nj, k % nk);
        return blocks(b)(c);
    }), nd::uivec(bi * ni, bj * nj, bk * nk));
}




//=============================================================================
template<typename P>
auto tile_blocks(std::map<nd::uivec_t<3>, nd::array_t<P, 3>> blocks, nd::uivec_t<3> shape)
{
    return tile_blocks(nd::make_array([blocks] (auto bi) { return blocks.at(bi); }, shape));
}




/**
 * @brief      Generate a sequence of 27 nd::uivec_t<3> items, in row-major
 *             order, identifying the block indexes isotropically surrounding a
 *             block at the given base index. Negative indexes are wrapped to
 *             the end (periodic topology), which is why this function requires
 *             an extent parameter.
 *
 * @param[in]  base_index  The base index, lying at the center of the 3 x 3 x 3
 *                         block
 * @param[in]  extent      The size of the block index space
 *
 * @return     A sequence of nd::uivec_t
 */
inline auto neighbors_27(nd::uivec_t<3> base_index, nd::uivec_t<3> extent)
{
    return seq::adapt(nd::index_space(3, 3, 3))
    | seq::map(to_numeric_array)
    | seq::map([i0=to_numeric_array(base_index), ex=to_numeric_array(extent)] (auto di)
    {
        return (i0 + di + ex - 1) % ex;
    })
    | seq::map(to_uivec);
}




/**
 * @brief      Generate a sequence of 9 nd::uivec_t<3> items, in row-major
 *             order, identifying the block indexes on a 3 x 3 square, oriented
 *             on one of the 3d axes, and centered on the given base index.
 *             Negative indexes are wrapped to the end (periodic topology),
 *             which is why this function requires an extent parameter.
 *
 * @param[in]  base_index  The base index, lying at the center of the 3 x 3
 *                         square
 * @param[in]  extent      The size of the block index space
 * @param[in]  axis        The axis normal to the plane containing the 3 x 3
 *                         square
 *
 * @return     A sequence of nd::uivec_t
 */
inline auto neighbors_9(nd::uivec_t<3> base_index, nd::uivec_t<3> extent, axis_3d axis)
{
    auto space = nd::index_space(axis == axis_3d::i ? 1 : 3, axis == axis_3d::j ? 1 : 3, axis == axis_3d::k ? 1 : 3);
    auto offset = numeric::array(axis == axis_3d::i ? 0 : 1, axis == axis_3d::j ? 0 : 1, axis == axis_3d::k ? 0 : 1);

    return seq::adapt(space)
    | seq::map(to_numeric_array)
    | seq::map([i0=to_numeric_array(base_index), ex=to_numeric_array(extent), offset] (auto di)
    {
        return (i0 + di + ex - offset) % ex;
    })
    | seq::map(to_uivec);
}




//=============================================================================
auto extend_periodic = [] (nd::uint count)
{
    return [count] (auto p)
    {
        return p | nd::extend_periodic(0, count) | nd::extend_periodic(1, count) | nd::extend_periodic(2, count);
    };
};

auto remove_surface = [] (nd::uint count)
{
    return [c=count] (auto x)
    {
        return nd::make_array(nd::indexing([x,c] (auto i, auto j, auto k)
        {
            return x(i + c, j + c, k + c);
        }), nd::uivec(shape(x, 0) - 2 * c, shape(x, 1) - 2 * c, shape(x, 2) - 2 * c));
    };
};

auto remove_transverse_i = [] (nd::uint count)
{
    return [c=count] (auto x)
    {
        return nd::make_array(nd::indexing([x, c] (auto i, auto j, auto k)
        {
            return x(i, j + c, k + c);
        }), nd::uivec(shape(x, 0), shape(x, 1) - 2 * c, shape(x, 2) - 2 * c));
    };
};

auto remove_transverse_j = [] (nd::uint count)
{
    return [c=count] (auto x)
    {
        return nd::make_array(nd::indexing([x, c] (auto i, auto j, auto k)
        {
            return x(i + c, j, k + c);
        }), nd::uivec(shape(x, 0) - 2 * c, shape(x, 1), shape(x, 2) - 2 * c));
    };
};

auto remove_transverse_k = [] (nd::uint count)
{
    return [c=count] (auto x)
    {
        return nd::make_array(nd::indexing([x, c] (auto i, auto j, auto k)
        {
            return x(i + c, j + c, k);
        }), nd::uivec(shape(x, 0) - 2 * c, shape(x, 1) - 2 * c, shape(x, 2)));
    };
};




//=============================================================================
template<typename P>
auto remove_transverse(nd::array_t<P, 3> array, unsigned count, axis_3d axis)
{
    return nd::make_array([array, count, axis] (auto i0)
    {
        return array(to_uivec(to_numeric_array(i0) + count * (1 - kronecker_delta(axis))));
    }, to_uivec(to_numeric_array(shape(array)) - 2 * count * (1 - kronecker_delta(axis))));
}

inline auto remove_transverse(unsigned count, axis_3d axis)
{
    return [count, axis] (auto array) { return remove_transverse(array, count, axis);};
}

} // namespace mesh




//=============================================================================
template<typename T>
auto to_string(const nd::shared_array<T, 3>& v)
{
    return "("
    + std::to_string(shape(v, 0)) + " "
    + std::to_string(shape(v, 1)) + " "
    + std::to_string(shape(v, 2)) + ")";
}

template<typename T>
auto to_string(const std::array<nd::shared_array<T, 3>, 3>& v)
{
    return to_string(v.at(0)) + " " + to_string(v.at(1)) + " " + to_string(v.at(2));
}




//=============================================================================
#ifdef DO_UNIT_TESTS
#include <iostream>
#include "core_unit_test.hpp"




//=============================================================================
inline void test_mesh_cartesian_3d()
{
    auto verts = mesh::unit_lattice(10, 10, 10);
    require(shape(verts) == nd::uivec(10, 10, 10));
    require(seq::to<std::vector>(mesh::neighbors_27(nd::uivec(2, 2, 2), nd::uivec(16, 16, 16))).at(0) == nd::uivec( 1,  1,  1));
    require(seq::to<std::vector>(mesh::neighbors_27(nd::uivec(2, 2, 2), nd::uivec(16, 16, 16))).at(5) == nd::uivec( 1,  2,  3));
    require(seq::to<std::vector>(mesh::neighbors_27(nd::uivec(0, 0, 0), nd::uivec(16, 16, 16))).at(0) == nd::uivec(15, 15, 15));
    require(seq::to<std::vector>(mesh::neighbors_27(nd::uivec(0, 0, 0), nd::uivec(16, 16, 16))).at(5) == nd::uivec(15,  0,  1));
}

#endif // DO_UNIT_TESTS
