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
#include <array>
#include <functional>




//=============================================================================
namespace numeric {

namespace detail {

template<typename T, typename U, std::size_t S, std::size_t... Is>
auto zip_arrays_impl(std::array<T, S> t, std::array<U, S> u, std::index_sequence<Is...>)
{
    return std::array{std::pair(std::get<Is>(t), std::get<Is>(u))...};
}

template<typename T, typename U, std::size_t S>
auto zip_arrays(std::array<T, S> t, std::array<U, S> u)
{
    return zip_arrays_impl(t, u, std::make_index_sequence<S>());
}

template<typename T, std::size_t S, typename FunctionType>
auto map_array(std::array<T, S> t, FunctionType fn)
{
    return std::apply([fn] (auto... ts) { return std::array{fn(ts)...}; }, t);
}

}




//=============================================================================
template<typename T, std::size_t S>
struct array_t
{
    const T& operator[](std::size_t i) const { return impl[i]; }
    T& operator[](std::size_t i) { return impl[i]; }
    const T& at(std::size_t i) const
    {
        if (i >= S)
            throw std::out_of_range("numeric::array_t (index out of range)");
        return impl[i];
    }
    T& at(std::size_t i)
    {
        if (i >= S)
            throw std::out_of_range("numeric::array_t (index out of range)");
        return impl[i];
    }
    std::array<T, S> impl;
};




//=============================================================================
template<typename... Ts>
auto array(Ts... ts)
{
    auto v = std::array{ts...};
    return array_t<typename decltype(v)::value_type, v.size()>{v};
}

template<typename T, std::size_t S>
auto array(std::array<T, S> t)
{
    return array_t<T, S>{t};
}

template<std::size_t S, typename T>
auto uniform(T value)
{
    auto v = array_t<T, S>{};

    for (std::size_t i = 0; i < S; ++i)
        v[i] = value;
    return v;
}

template<std::size_t S, typename T=unsigned long>
auto range()
{
    auto v = array_t<T, S>{};

    for (std::size_t i = 0; i < S; ++i)
        v[i] = i;
    return v;
}

template<std::size_t I, typename T, std::size_t S, typename = std::enable_if_t<I < S>>
const auto& get(const array_t<T, S>& t)
{
    return std::get<I>(t.impl);
}

template<std::size_t I, typename T, std::size_t S, typename = std::enable_if_t<I < S>>
auto& get(array_t<T, S>& t)
{
    return std::get<I>(t.impl);
}

template<typename T, std::size_t S>
auto as_tuple(array_t<T, S> a)
{
    return a.impl;
}

template<typename T, std::size_t S, typename FunctionType>
auto map(array_t<T, S> t, FunctionType f)
{   
    return array(detail::map_array(t.impl, f));
}

template<typename T, typename U, std::size_t S>
auto zip(array_t<T, S> t, array_t<U, S> u)
{
    return array(detail::zip_arrays(t.impl, u.impl));
}

template<typename T, typename U, std::size_t S>
auto apply_to(array_t<T, S> t, array_t<U, S> u)
{
    return map(zip(t, u), [] (auto tu) { return get<0>(tu)(get<1>(tu)); });
}

template<typename T, std::size_t S>
bool any(array_t<T, S> t)
{
    for (auto a : t.impl)
        if (a)
            return true;
    return false;
}

template<typename T, std::size_t S>
bool all(array_t<T, S> t)
{
    for (auto a : t.impl)
        if (! a)
            return false;
    return true;
}

template<typename T, std::size_t S, typename FunctionType>
auto reduce(array_t<T, S> t, T seed, FunctionType f)
{   
    for (std::size_t i = 0; i < S; ++i)
        seed = f(seed, t[i]);
    return seed;
}

template<typename T, std::size_t S>
auto sum(array_t<T, S> t)
{   
    return reduce(t, T(0), std::plus<>());
}

template<typename T, std::size_t S>
auto product(array_t<T, S> t)
{   
    return reduce(t, T(1), std::multiplies<>());
}

template<typename T, std::size_t S>
auto min(array_t<T, S> t)
{   
    return reduce(t, t.at(0), [] (auto a, auto b) { return std::min(a, b); });
}

template<typename T, std::size_t S>
auto max(array_t<T, S> t)
{   
    return reduce(t, t.at(0), [] (auto a, auto b) { return std::max(a, b); });
}

template<typename T, std::size_t S>
auto replace(array_t<T, S> t, std::size_t index, T value)
{   
    t.at(index) = value;
    return t;
}

template<typename T, std::size_t S, typename FunctionType>
auto update(array_t<T, S> t, std::size_t index, FunctionType f)
{
    return replace(t, index, f(t.at(index)));
}

template<typename T, typename U, std::size_t S>
auto construct(numeric::array_t<U, S> a)
{
    return map(a, [] (auto x) { return T(x); });
}




/**
 * @brief      Return a boolean sequence {a} representing a number:
 *
 *             value = a[0] * 2^0 + a[1] * 2^1 + ...
 *
 * @param[in]  value     The number to represent
 *
 * @tparam     BitCount  The number of bits to include
 *
 * @return     A boolean sequence, with the most significant bit at the
 *             beginning.
 */
template<unsigned long BitCount>
auto binary_repr(std::size_t value)
{
    return map(numeric::range<BitCount>(), [value] (auto n) { return bool(value & (1 << n)); });
}




/**
 * @brief      Turn a binary representation of a number into a 64-bit unsigned
 *             integer.
 *
 * @param[in]  bits      The bits (as returned from binary_repr)
 *
 * @tparam     BitCount  The number of bits
 *
 * @return     The decimal representation
 */
template<unsigned long BitCount>
unsigned long to_integral(numeric::array_t<bool, BitCount> bits)
{
    return sum(apply_to(map(numeric::range<BitCount>(), [] (auto e) { return [e] (bool y) { return (1 << e) * y; }; }), bits));
}




//=============================================================================
template<typename T, typename U, std::size_t S> auto operator+ (array_t<T, S> a, array_t<U, S> b) { return map(zip(a, b), [] (auto t) { return std::apply(std::plus<>(), t); }); }
template<typename T, typename U, std::size_t S> auto operator- (array_t<T, S> a, array_t<U, S> b) { return map(zip(a, b), [] (auto t) { return std::apply(std::minus<>(), t); }); }
template<typename T, typename U, std::size_t S> auto operator* (array_t<T, S> a, array_t<U, S> b) { return map(zip(a, b), [] (auto t) { return std::apply(std::multiplies<>(), t); }); }
template<typename T, typename U, std::size_t S> auto operator/ (array_t<T, S> a, array_t<U, S> b) { return map(zip(a, b), [] (auto t) { return std::apply(std::divides<>(), t); }); }
template<typename T, typename U, std::size_t S> auto operator% (array_t<T, S> a, array_t<U, S> b) { return map(zip(a, b), [] (auto t) { return std::apply(std::modulus<>(), t); }); }
template<typename T, typename U, std::size_t S> auto operator==(array_t<T, S> a, array_t<U, S> b) { return all(map(zip(a, b), [] (auto t) { return std::apply(std::equal_to<>(), t); })); }
template<typename T, typename U, std::size_t S> auto operator!=(array_t<T, S> a, array_t<U, S> b) { return any(map(zip(a, b), [] (auto t) { return std::apply(std::not_equal_to<>(), t); })); }

template<typename T, typename U, std::size_t S> auto operator+(array_t<T, S> a, U b) { return a + uniform<S>(b); }
template<typename T, typename U, std::size_t S> auto operator-(array_t<T, S> a, U b) { return a - uniform<S>(b); }
template<typename T, typename U, std::size_t S> auto operator*(array_t<T, S> a, U b) { return a * uniform<S>(b); }
template<typename T, typename U, std::size_t S> auto operator/(array_t<T, S> a, U b) { return a / uniform<S>(b); }
template<typename T, typename U, std::size_t S> auto operator%(array_t<T, S> a, U b) { return a % uniform<S>(b); }
template<typename T, typename U, std::size_t S> auto operator+(U b, array_t<T, S> a) { return uniform<S>(b) + a; }
template<typename T, typename U, std::size_t S> auto operator-(U b, array_t<T, S> a) { return uniform<S>(b) - a; }
template<typename T, typename U, std::size_t S> auto operator*(U b, array_t<T, S> a) { return uniform<S>(b) * a; }
template<typename T, typename U, std::size_t S> auto operator/(U b, array_t<T, S> a) { return uniform<S>(b) / a; }
template<typename T, typename U, std::size_t S> auto operator%(U b, array_t<T, S> a) { return uniform<S>(b) % a; }

template<typename T, std::size_t S> auto operator+(array_t<T, S> a) { return map(a, [] (auto t) { return +t; }); }
template<typename T, std::size_t S> auto operator-(array_t<T, S> a) { return map(a, [] (auto t) { return -t; }); }

} // namespace numeric




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_numeric_array()
{
    require(numeric::array(1, 2) == numeric::array(1.0, 2.0));
    require(numeric::array(1, 2) * 2.1 == numeric::array(2.1, 4.2));
    require(sum(numeric::array(1, 2, 3, 4)) == 10);
    require(product(numeric::array(1., 2., 3., 4.)) == 24);

    require(to_integral(numeric::array(false, false, false)) == 0);
    require(to_integral(numeric::array(true, true, true)) == 7);
    require(to_integral(numeric::array(true, false, false)) == 1); // the 2^0 bit is on the left
    require(to_integral(numeric::binary_repr<3>(6)) == 6);
}

#endif // DO_UNIT_TESTS
