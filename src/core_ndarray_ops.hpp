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
#include "core_ndarray.hpp"




//=============================================================================
namespace nd {




//=============================================================================
template<typename ArrayType>
auto multiply(ArrayType b)
{
    return [b] (auto a)
    {
        return a * b;
    };
}

template<typename ArrayType>
auto divide(ArrayType b)
{
    return [b] (auto a)
    {
        return a / b;
    };
}

inline auto adjacent_zip(uint axis=0)
{
    return [axis] (auto x)
    {
        return zip(select(x, axis, 0, -1), select(x, axis, 1));
    };
}

inline auto adjacent_zip3(uint axis=0)
{
    return [axis] (auto x)
    {
        return zip(select(x, axis, 0, -2), select(x, axis, 1, -1), select(x, axis, 2));
    };
}

inline auto adjacent_zip4(uint axis=0)
{
    return [axis] (auto x)
    {
        return zip(select(x, axis, 0, -3), select(x, axis, 1, -2), select(x, axis, 2, -1), select(x, axis, 3));
    };
}

inline auto adjacent_mean(uint axis=0)
{
    return [axis] (auto x)
    {
        return 0.5 * (select(x, axis, 0, -1) + select(x, axis, 1));
    };
}

inline auto adjacent_diff(uint axis=0)
{
    return [axis] (auto x)
    {
        return select(x, axis, 1) - select(x, axis, 0, -1);
    };
}

inline auto extend_periodic(uint axis=0, uint count=1)
{
    return [count, axis] (auto x)
    {
        return select(x, axis, -count)
        | nd::concat(x, axis)
        | nd::concat(select(x, axis, 0, count), axis);
    };
}

inline auto extend_zero_gradient_lower(uint axis=0, uint count=1)
{
    return [axis, count] (auto x)
    {
        return nd::concat(repeat(select(x, axis, 0, 1), axis, count), x, axis);
    };
}

inline auto extend_zero_gradient_upper(uint axis=0, uint count=1)
{
    return [axis, count] (auto x)
    {
        return nd::concat(x, repeat(select(x, axis, -1), axis, count), axis);
    };
}

inline auto extend_zero_gradient(uint axis=0, uint count=1)
{
    return [count, axis] (auto x)
    {
        return x | extend_zero_gradient_lower(axis, count) | extend_zero_gradient_upper(axis, count);
    };
}

template<typename T>
inline auto extend_uniform_lower(T value, uint axis=0, uint count=1)
{
    return [count, axis, value] (auto x)
    {
        return nd::concat(nd::uniform(value, replace(shape(x), axis, count)), x, axis);
    };
}

template<typename T>
auto extend_uniform_upper(T value, uint axis=0, uint count=1)
{
    return [count, axis, value] (auto x)
    {
        return nd::concat(x, nd::uniform(value, replace(shape(x), axis, count)), axis);
    };
}

template<typename T>
auto extend_uniform(T value, uint axis=0, uint count=1)
{
    return [count, axis, value] (auto x)
    {
        return x | extend_uniform_lower(value, axis, count) | extend_uniform_upper(value, axis, count);
    };
}

inline auto extend_zeros(uint axis=0, uint count=1)
{
    return [count, axis] (auto x)
    {
        return x | extend_uniform(typename decltype(x)::value_type(), axis, count);
    };
}

template<typename T>
auto construct()
{
    return [] (auto a)
    {
        return map(a, [] (auto x) { return T(x); });
    };
}




/**
 * @brief      Remove a partition from a pair of bin edges and corresponding
 *             masses. This operation removes the edge at index i, and combines
 *             the masses in bins i - 1 and i. The size of the edge array must
 *             be one larger than that of the mass array, and the first and last
 *             edges cannot be removed.
 *
 * @param[in]  edge_array     The array of bin edges
 * @param[in]  mass_array     The array of bin masses
 * @param[in]  edge_index     The edge index to remove
 *
 * @tparam     ProviderType1  The type of the edge_array provider
 * @tparam     ProviderType2  The type of the mass_array provider
 *
 * @return     A pair (new_edge_array, new_mass_array), each shortened by a
 *             single element.
 */
template<typename ProviderType1, typename ProviderType2>
auto remove_partition(array_t<ProviderType1, 1> edge_array, array_t<ProviderType2, 1> mass_array, uint edge_index)
{
    if (size(edge_array) != size(mass_array) + 1)
        throw std::invalid_argument("nd::remove_partition (size(edge_array) != size(mass_array) + 1)");

    if (edge_index == 0 || edge_index >= size(edge_array) - 1)
        throw std::out_of_range("nd::remove_partition (edge_index out of range)");

    auto new_edge_array = select(edge_array, 0, 0, edge_index)
    | concat(select(edge_array, 0, edge_index + 1));

    auto new_mass_array = select(mass_array, 0, 0, edge_index - 1)
    | concat(nd::from(mass_array(edge_index - 1) + mass_array(edge_index)))
    | concat(select(mass_array, 0, edge_index + 1));

    return std::pair(new_edge_array, new_mass_array);
}

} // namespace nd




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_ndarray_ops()
{
    {
        auto [e, m] = nd::remove_partition(nd::range(4), nd::from(1.0, 2.0, 3.0), 1);
        require(all(e == nd::from(0, 2, 3)));
        require(all(m == nd::from(3.0, 3.0)));
    }

    {
        auto [e, m] = nd::remove_partition(nd::range(4), nd::from(1.0, 2.0, 3.0), 2);
        require(all(e == nd::from(0, 1, 3)));
        require(all(m == nd::from(1.0, 5.0)));
    }
}

#endif // DO_UNIT_TESTS
