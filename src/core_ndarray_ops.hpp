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

inline auto extend_extrap()
{
    return [] (auto x)
    {
        auto dx_l = x(1) - x(0);
        auto dx_r = x(size(x) - 1) - x(size(x) - 2);
        auto xl = x(0)           - dx_l;
        auto xr = x(size(x) - 1) + dx_r;
        return from(xl) | concat(x) | concat(from(xr));
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
 * @brief      Convert a 1d array of 1d arrays of T into a 1d array of T, by
 *             joining the arrays end-to-end.
 *
 * @param[in]  a          The array of arrays
 *
 * @tparam     P          The provider type of a
 * @tparam     T          The array type held by a
 * @tparam     <unnamed>  enable-if T is an array
 * @tparam     <unnamed>  enable-if T has rank 1
 *
 * @return     A 1d shared array of value type T
 */
template<
    typename P,
    typename T = typename nd::array_t<P, 1>::value_type,
    typename = std::enable_if_t<nd::is_array<T>::value>,
    typename = std::enable_if_t<nd::array_t<P, 1>::rank == 1>>
auto flat(nd::array_t<P, 1> a)
{
    auto total_cells = map(a, [] (auto ai) { return size(ai); }) | nd::sum();
    auto result = nd::make_unique_array<typename T::value_type>(nd::uivec(total_cells));
    auto n = nd::uint(0);

    for (std::size_t i = 0; i < size(a); ++i)
        for (std::size_t j = 0; j < size(a(i)); ++j)
            result(n++) = a(i)(j);

    return nd::make_shared_array(std::move(result));
}

inline auto flat()
{
    return [] (auto a) { return flat(a); };
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

    if (edge_index == 0 || edge_index + 1 >= size(edge_array))
        throw std::out_of_range("nd::remove_partition (edge_index out of range)");

    auto new_edge_array = select(edge_array, 0, 0, edge_index)
    | concat(select(edge_array, 0, edge_index + 1))
    | nd::to_shared();

    auto new_mass_array = select(mass_array, 0, 0, edge_index - 1)
    | concat(nd::from(mass_array(edge_index - 1) + mass_array(edge_index)))
    | concat(select(mass_array, 0, edge_index + 1))
    | nd::to_shared();

    return std::pair(new_edge_array, new_mass_array);
}




/**
 * @brief      Add a new partition, bisecting the bin at the given index and
 *             distributing its mass equally.
 *
 * @param[in]  edge_array     The array of bin edges
 * @param[in]  mass_array     The array of bin masses
 * @param[in]  bin_index      The bin to subdivide
 *
 * @tparam     ProviderType1  The type of the edge_array provider
 * @tparam     ProviderType2  The type of the mass_array provider
 *
 * @return     A pair (new_edge_array, new_mass_array), each lengthened by a
 *             single element.
 */
template<typename ProviderType1, typename ProviderType2>
auto add_partition(array_t<ProviderType1, 1> edge_array, array_t<ProviderType2, 1> mass_array, uint bin_index)
{
    if (size(edge_array) != size(mass_array) + 1)
        throw std::invalid_argument("nd::add_partition (size(edge_array) != size(mass_array) + 1)");

    if (bin_index + 1 > size(mass_array))
        throw std::out_of_range("nd::add_partition (bin_index out of range)");

    auto new_edge_array = select(edge_array, 0, 0, bin_index + 1)
    | concat(nd::from(0.5 * (edge_array(bin_index) + edge_array(bin_index + 1))))
    | concat(select(edge_array, 0, bin_index + 1))
    | nd::to_shared();

    auto new_mass_array = select(mass_array, 0, 0, bin_index)
    | concat(nd::from(0.5 * mass_array(bin_index), 0.5 * mass_array(bin_index)))
    | concat(select(mass_array, 0, bin_index + 1))
    | nd::to_shared();

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

    {
        auto [e, m] = nd::add_partition(nd::linspace(0.0, 1.0, 5), nd::from(1.0, 1.0, 1.0, 1.0), 0);
        require(all(e == nd::from(0.0, 0.125, 0.25, 0.50, 0.75, 1.0)));
        require(all(m == nd::from(0.5, 0.5, 1.0, 1.0, 1.0)));
    }

    {
        auto [e, m] = nd::add_partition(nd::linspace(0.0, 1.0, 5), nd::from(1.0, 1.0, 1.0, 1.0), 1);
        require(all(e == nd::from(0.0, 0.25, 0.375, 0.50, 0.75, 1.0)));
        require(all(m == nd::from(1.0, 0.5, 0.5, 1.0, 1.0)));
    }
}

#endif // DO_UNIT_TESTS
