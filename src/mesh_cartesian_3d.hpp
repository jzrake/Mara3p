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
#include <vector>
#include "core_geometric.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_numeric_array.hpp"
#include "core_util.hpp"




//=============================================================================
namespace mesh {




/**
 * @brief      This enum is used to represent an axis in a 3d array
 */
enum class axis_3d {
    i, j, k
};

using block_index_3d_t = numeric::array_t<unsigned long, 3>;




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
template<typename P>
auto tile_blocks_27(std::vector<nd::array_t<P, 3>> blocks)
{
    if (blocks.size() != 27)
        throw std::invalid_argument("mesh::tile_blocks_27 (must supply exactly 27 blocks)");

    for (const auto& block : blocks)
        if (shape(block) != shape(blocks.at(0)))
            throw std::invalid_argument("mesh::tile_blocks_27 (supplied blocks have non-uniform shapes)");

    auto ni = shape(blocks.at(0), 0);
    auto nj = shape(blocks.at(0), 1);
    auto nk = shape(blocks.at(0), 2);

    return nd::make_array(nd::indexing([ni, nj, nk, blocks] (auto i, auto j, auto k)
    {
        auto s = nd::uivec(1, 3, 9);
        auto b = nd::uivec(i / ni, j / nj, k / nk);
        auto c = nd::uivec(i % ni, j % nj, k % nk);
        return blocks[dot(s, b)](c);
    }), nd::uivec(3 * ni, 3 * nj, 3 * nk));
}




//=============================================================================
template<typename P>
auto tile_blocks_9(std::vector<nd::array_t<P, 3>> blocks, axis_3d axis)
{
    if (blocks.size() != 9)
        throw std::invalid_argument("mesh::tile_blocks_9 (must supply exactly 9 blocks)");

    for (const auto& block : blocks)
        if (shape(block) != shape(blocks.at(0)))
            throw std::invalid_argument("mesh::tile_blocks_9 (supplied blocks have non-uniform shapes)");

    auto ni = shape(blocks.at(0), 0);
    auto nj = shape(blocks.at(0), 1);
    auto nk = shape(blocks.at(0), 2);

    return nd::make_array(nd::indexing([ni, nj, nk, blocks, axis] (auto i, auto j, auto k)
    {
        switch (axis)
        {
            case axis_3d::i: {
                auto s = nd::uivec(1, 3);
                auto b = nd::uivec(j / nj, k / nk);
                auto c = nd::uivec(i, j % nj, k % nk);
                return blocks[dot(s, b)](c);                
            }
            case axis_3d::j: {
                auto s = nd::uivec(1, 3);
                auto b = nd::uivec(k / nk, i / ni);
                auto c = nd::uivec(i % ni, j, k % nk);
                return blocks[dot(s, b)](c);                
            }
            case axis_3d::k: {
                auto s = nd::uivec(1, 3);
                auto b = nd::uivec(i / ni, j / nj);
                auto c = nd::uivec(i % ni, j % nj, k);
                return blocks[dot(s, b)](c);                
            }
        }
    }), nd::uivec(
        (axis == axis_3d::i ? 1 : 3) * ni,
        (axis == axis_3d::j ? 1 : 3) * nj,
        (axis == axis_3d::k ? 1 : 3) * nk));
}

template<typename P> auto tile_blocks_faces_9i(std::vector<nd::array_t<P, 3>> b) { return tile_blocks_9(b, axis_3d::i); }
template<typename P> auto tile_blocks_faces_9j(std::vector<nd::array_t<P, 3>> b) { return tile_blocks_9(b, axis_3d::j); }
template<typename P> auto tile_blocks_faces_9k(std::vector<nd::array_t<P, 3>> b) { return tile_blocks_9(b, axis_3d::k); }




//=============================================================================
inline std::vector<block_index_3d_t> neighbors_27(block_index_3d_t index, block_index_3d_t extent)
{
    auto wi = [e=extent[0]] (auto i) { return (i + e) % e; };
    auto wj = [e=extent[1]] (auto j) { return (j + e) % e; };
    auto wk = [e=extent[2]] (auto k) { return (k + e) % e; };
    auto stride = nd::uivec(1, 3, 9);
    auto result = std::vector<block_index_3d_t>(27);

    for (auto d : nd::index_space(nd::uivec(3, 3, 3)))
    {
        result[dot(stride, d)] = numeric::array(
            wi(index[0] + d[0] - 1),
            wj(index[1] + d[1] - 1),
            wk(index[2] + d[2] - 1));
    }
    return result;
}




//=============================================================================
inline std::vector<block_index_3d_t> neighbors_9(block_index_3d_t index, block_index_3d_t extent, axis_3d axis)
{
    auto wi = [e=extent[0]] (auto i) { return (i + e) % e; };
    auto wj = [e=extent[1]] (auto j) { return (j + e) % e; };
    auto wk = [e=extent[2]] (auto k) { return (k + e) % e; };
    auto [i, j, k] = index.impl;
    auto stride = nd::uivec(1, 3);
    auto result = std::vector<block_index_3d_t>(9);

    for (auto d : nd::index_space(nd::uivec(3, 3)))
    {
        auto [da, db] = d.impl;

        switch (axis)
        {
            case axis_3d::i: result[dot(stride, d)] = numeric::array(i, wj(j + da - 1), wk(k + db - 1)); break;
            case axis_3d::j: result[dot(stride, d)] = numeric::array(wi(i + db - 1), j, wk(k + da - 1)); break;
            case axis_3d::k: result[dot(stride, d)] = numeric::array(wi(i + da - 1), wk(k + db - 1), k); break;
        }
    }
    return result;
}




//=============================================================================
auto extend_periodic = [] (nd::uint count)
{
    return [count] (auto p)
    {
        return p | nd::extend_periodic(0, count) | nd::extend_periodic(1, count) | nd::extend_periodic(2, count);
    };
};




//=============================================================================
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
        return nd::make_array(nd::indexing([x,c] (auto i, auto j, auto k)
        {
            return x(i, j + c, k + c);
        }), nd::uivec(shape(x, 0), shape(x, 1) - 2 * c, shape(x, 2) - 2 * c));
    };
};

auto remove_transverse_j = [] (nd::uint count)
{
    return [c=count] (auto x)
    {
        return nd::make_array(nd::indexing([x,c] (auto i, auto j, auto k)
        {
            return x(i + c, j, k + c);
        }), nd::uivec(shape(x, 0) - 2 * c, shape(x, 1), shape(x, 2) - 2 * c));
    };
};

auto remove_transverse_k = [] (nd::uint count)
{
    return [c=count] (auto x)
    {
        return nd::make_array(nd::indexing([x,c] (auto i, auto j, auto k)
        {
            return x(i + c, j + c, k);
        }), nd::uivec(shape(x, 0) - 2 * c, shape(x, 1) - 2 * c, shape(x, 2)));
    };
};

} // namespace mesh




//=============================================================================
#ifdef DO_UNIT_TESTS
#include <iostream>
#include "core_unit_test.hpp"




//=============================================================================
inline void test_mesh_cartesian_3d()
{
    auto verts = mesh::unit_lattice(10, 10, 10);
    require(shape(verts) == nd::uivec(10, 10, 10));

    {
        auto blocks = std::vector<nd::shared_array<geometric::euclidean_vector_t<double>, 3>>();
        auto A = mesh::unit_lattice(16, 17, 18) | nd::to_shared();

        for (int i = 0; i < 27; ++i)
            blocks.push_back(A);

        auto B = mesh::tile_blocks_27(blocks);
        require(shape(B) == nd::uivec(16 * 3, 17 * 3, 18 * 3));
        require(B(0, 0, 0) == B(16, 17, 18));
        require(B(16, 17, 18) == B(32, 34, 36));
        require(B(16, 17, 19) == B(32, 34, 37));
    }

    require(mesh::neighbors_27(numeric::array(2UL, 2UL, 2UL), numeric::array(16UL, 16UL, 16UL)).at(0) == numeric::array( 1UL,  1UL, 1UL));
    require(mesh::neighbors_27(numeric::array(2UL, 2UL, 2UL), numeric::array(16UL, 16UL, 16UL)).at(5) == numeric::array( 3UL,  2UL, 1UL));
    require(mesh::neighbors_27(numeric::array(0UL, 0UL, 0UL), numeric::array(16UL, 16UL, 16UL)).at(0) == numeric::array(15UL, 15UL, 15UL));
    require(mesh::neighbors_27(numeric::array(0UL, 0UL, 0UL), numeric::array(16UL, 16UL, 16UL)).at(5) == numeric::array( 1UL,  0UL, 15UL));
}

#endif // DO_UNIT_TESTS
