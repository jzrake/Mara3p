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
#include <tuple>
#include <functional>




//=============================================================================
namespace numeric {

namespace detail {

template<typename... Ts, typename... Us, std::size_t... Is>
auto zip_tuples_impl(std::tuple<Ts...> t, std::tuple<Us...> u, std::index_sequence<Is...>)
{
    return std::tuple(std::pair(std::get<Is>(t), std::get<Is>(u))...);
}

template<typename... Ts, typename... Us>
auto zip_tuples(std::tuple<Ts...> t, std::tuple<Us...> u)
{
    return zip_tuples_impl(t, u, std::make_index_sequence<sizeof...(Ts)>());
}

template<typename FunctionType, typename... Ts>
auto map_tuple(std::tuple<Ts...> t, FunctionType fn)
{
    return std::apply([fn] (auto... ts) { return std::tuple(fn(ts)...); }, t);
}

template<typename FunctionType>
auto apply_of(FunctionType f)
{
    return [f] (auto t) { return std::apply(f, t); };
}

} // namespace detail




//=============================================================================
template<typename... Ts>
struct tuple_t
{
    tuple_t(std::tuple<Ts...> impl) : impl(impl) {}
    tuple_t() {}
    std::tuple<Ts...> impl;
};




//=============================================================================
template<typename... Ts>
auto tuple(Ts... ts)
{
    return tuple_t<Ts...>{{ts...}};
}

template<typename... Ts>
auto tuple(std::tuple<Ts...> t)
{
    return tuple_t<Ts...>{t};
}

template<typename... TupleTypes>
auto tuple_cat(TupleTypes... ts)
{
    return tuple(std::tuple_cat(ts.impl...));
}




//=============================================================================
template<std::size_t I, typename... Ts, typename = std::enable_if_t<I < sizeof...(Ts)>>
const auto& get(const tuple_t<Ts...>& t)
{
    return std::get<I>(t.impl);
}

template<std::size_t I, typename... Ts, typename = std::enable_if_t<I < sizeof...(Ts)>>
auto& get(tuple_t<Ts...>& t)
{
    return std::get<I>(t.impl);
}

template<typename... Ts, typename FunctionType>
auto map(tuple_t<Ts...> t, FunctionType f)
{
    return tuple(detail::map_tuple(t.impl, f));
}

template<typename... Ts, typename... Us, typename = std::enable_if_t<sizeof...(Ts) == sizeof...(Us)>>
auto zip(tuple_t<Ts...> a, tuple_t<Us...> b)
{
    return tuple(detail::zip_tuples(a.impl, b.impl));
}

template<typename... Ts>
bool any(tuple_t<Ts...> t)
{
    auto f = [] (auto... as)
    {
        for (auto a : {bool(as)...})
            if (a)
                return true;
        return false;
    };
    return std::apply(f, t.impl);    
}

template<typename... Ts>
bool all(tuple_t<Ts...> t)
{
    auto f = [] (auto... as)
    {
        for (auto a : {bool(as)...})
            if (! a)
                return false;
        return true;
    };
    return std::apply(f, t.impl);    
}




//=============================================================================
template<typename... Ts, typename... Us> auto operator+(tuple_t<Ts...> a, tuple_t<Us...> b) { return tuple(detail::map_tuple(detail::zip_tuples(a.impl, b.impl), detail::apply_of(std::plus<>()))); }
template<typename... Ts, typename... Us> auto operator-(tuple_t<Ts...> a, tuple_t<Us...> b) { return tuple(detail::map_tuple(detail::zip_tuples(a.impl, b.impl), detail::apply_of(std::minus<>()))); }
template<typename... Ts, typename U> auto operator*(tuple_t<Ts...> a, U b) { return tuple(detail::map_tuple(a.impl, [b] (auto ai) { return ai * b; })); }
template<typename... Ts, typename U> auto operator/(tuple_t<Ts...> a, U b) { return tuple(detail::map_tuple(a.impl, [b] (auto ai) { return ai / b; })); }
template<typename... Ts, typename U> auto operator*(U b, tuple_t<Ts...> a) { return tuple(detail::map_tuple(a.impl, [b] (auto ai) { return b * ai; })); }
template<typename... Ts, typename U> auto operator/(U b, tuple_t<Ts...> a) { return tuple(detail::map_tuple(a.impl, [b] (auto ai) { return b / ai; })); }
template<typename... Ts, typename... Us> auto operator==(tuple_t<Ts...> a, tuple_t<Us...> b) { return all(tuple(detail::map_tuple(detail::zip_tuples(a.impl, b.impl), detail::apply_of(std::equal_to<>())))); }
template<typename... Ts, typename... Us> auto operator!=(tuple_t<Ts...> a, tuple_t<Us...> b) { return any(tuple(detail::map_tuple(detail::zip_tuples(a.impl, b.impl), detail::apply_of(std::not_equal_to<>())))); }
template<typename... Ts> auto operator+(tuple_t<Ts...> a) { return map(a, [] (auto t) { return +t; }); }
template<typename... Ts> auto operator-(tuple_t<Ts...> a) { return map(a, [] (auto t) { return -t; }); }

} // namespace numeric




//=============================================================================
#ifdef DO_UNIT_TESTS
#include <cmath>
#include "core_unit_test.hpp"




//=============================================================================
inline void test_numeric_tuple()
{
    auto a = numeric::tuple(1, 2.3, 3.f);
    auto b = map(a, [] (auto x) { return std::sqrt(x); });

    static_assert(std::is_same_v<std::decay_t<decltype(numeric::get<0>(b))>, double>);
    static_assert(std::is_same_v<std::decay_t<decltype(numeric::get<1>(b))>, double>);
    static_assert(std::is_same_v<std::decay_t<decltype(numeric::get<2>(b))>, float>);

    require(numeric::tuple(1, 2) == numeric::tuple(1.0, 2.0));
    require(numeric::tuple(1, 2) != numeric::tuple(1.0, 2.1));
    require(numeric::tuple(1, 2) * 2.1 == numeric::tuple(2.1, 4.2));
}

#endif // DO_UNIT_TESTS
