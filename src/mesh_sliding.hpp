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
#include <optional>
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"




//=============================================================================
namespace mesh {


//=============================================================================
template<typename P, typename Q>
auto bin_values(nd::array_t<P, 1> masses, nd::array_t<Q, 1> bin_indexes, unsigned num_bins)
{
    auto result = nd::make_unique_array<typename nd::array_t<P, 1>::value_type>(nd::uivec(num_bins));

    for (auto [mass, bin] : zip(masses, bin_indexes))
    {
        if (bin >= num_bins)
        {
            continue;
            // throw std::out_of_range("mesh::bin_values");
        }
        result(bin) = result(bin) + mass;
    }
    return nd::make_shared_array(std::move(result));
}




//=============================================================================
template<typename P>
auto merge_sort(nd::array_t<P, 1> L, nd::array_t<P, 1> R)
{
    using value_type = typename nd::array_t<P, 1>::value_type;

    auto result = nd::make_unique_array<value_type>(nd::uivec(size(L) + size(R)));
    auto il     = nd::uint(0);
    auto ir     = nd::uint(0);
    auto n      = nd::uint(0);

    while (il < size(L) || ir < size(R))
    {
        if (ir == size(R) || L(il) <= R(ir))
        {
            result(n++) = L(il++);
        }
        if (il == size(L) || R(ir) <= L(il))
        {
            result(n++) = R(ir++);
        }
    }
    return nd::make_shared_array(std::move(result));
}




//=============================================================================
template<typename PositionType>
struct face_description_t
{
    std::optional<nd::uint> il; // index of zone on the left
    std::optional<nd::uint> ir; // index of zone on the right
    PositionType leading;  // coordinate of the vertex at the face's leading edge
    PositionType trailing; // coordinate of the vertex at the face's trailing edge
};




/**
 * @brief      Construct segments lying along the surface dividing two rows of
 *             quadrilateral cells.
 *
 * @param[in]  L     The positions of the longitudinal faces on the left row
 * @param[in]  R     The positions of the longitudinal faces on the right row
 *
 * @tparam     P     The provider type of the input arrays
 *
 * @return     An array of faces
 *
 * @note       This function assumes that the positions in both arrays are
 *             monotonically increasing. It also throws an exception if either
 *             row has a cell that does not overlap the neighbor.
 */
template<typename P>
auto transverse_faces(nd::array_t<P, 1> L, nd::array_t<P, 1> R)
{
    /*                 
     *
     *             |   x---|
     *             |   |   |  ]- transverse face [null / 3]
     *             |---x   |
     *             |   |   |
     *             |   x---|  <- longitudinal face on the right
     *             |   |   |
     *             |---x   |  <- longitudinal face on the left
     *             |   |   |
     *             |   x---|
     *             |   |   |
     *             |   x---|
     *             |   |   |
     *             |---x   |
     *             |   |   |  ]- transverse face [0 / 0]
     *             |   x---|
     *             |   |   |  ]- transverse face [0 / null]
     *             |---x   |
     */
    using face_type = face_description_t<typename nd::array_t<P, 1>::value_type>;

    auto result = nd::make_unique_array<face_type>(nd::uivec(size(L) + size(R) - 1));
    auto Lm = zip(L, nd::uniform('L', shape(L)));
    auto Rm = zip(R, nd::uniform('R', shape(R)));
    auto il = -1;
    auto ir = -1;
    auto n = nd::uint(0);

    auto check = [] (int i, nd::uint size) -> std::optional<nd::uint>
    {
        if (i < 0 || i + 1 >= size)
        {
            return std::nullopt;
        }
        return nd::uint(i);
    };

    for (auto [m0, m1] : merge_sort(Lm, Rm) | nd::adjacent_zip())
    {
        result(n++) = {
            check(il += (std::get<1>(m0) == 'L'), size(L)),
            check(ir += (std::get<1>(m0) == 'R'), size(R)),
            std::get<0>(m0),
            std::get<0>(m1),
        };
    }
    return nd::make_shared_array(std::move(result));
}

} // namespace mesh




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_mesh_sliding()
{
    auto L = nd::linspace(0.00, 1.00, 11);
    auto R = nd::linspace(0.05, 1.05, 11);
    auto F = mesh::transverse_faces(L, R);

    require(size(F) == 21);
    require(front(F).il == 0);
    require(front(F).ir == std::nullopt);
    require(back(F).il == std::nullopt);
    require(back(F).ir == 9);
}

#endif // DO_UNIT_TESTS