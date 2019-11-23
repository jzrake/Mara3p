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
#include "core_geometric.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"




//=============================================================================
namespace mesh {




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
template<typename T>
auto solenoidal_difference(nd::array_t<T, 3> e1, nd::array_t<T, 3> e2, nd::array_t<T, 3> e3)
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
 * @tparam     T     The provider type of the arrays
 *
 * @return     A single 3D array
 */
template<typename T>
auto divergence_difference(nd::array_t<T, 3> f1, nd::array_t<T, 3> f2, nd::array_t<T, 3> f3)
{
    auto d1 = nd::adjacent_diff(0);
    auto d2 = nd::adjacent_diff(1);
    auto d3 = nd::adjacent_diff(2);
    return d1(f1) + d2(f2) + d3(f3);
}

} // namespace mesh




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_mesh_cartesian_3d()
{
    auto verts = mesh::unit_lattice(10, 10, 10);
    require(shape(verts) == nd::uivec(10, 10, 10));
}

#endif // DO_UNIT_TESTS
