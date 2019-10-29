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




//=============================================================================
namespace numeric {




//=============================================================================
template<typename T>
struct optional_t
{
    std::optional<T> impl;
};




//=============================================================================
template<typename T>
auto optional(T t)
{
    return optional_t<T>{{t}};
}

template<typename T>
auto optional()
{
    return optional_t<T>{{}};
}

template<typename T>
bool has_value(const optional_t<T>& o)
{
    return o.impl.has_value();
}

template<typename T>
const T& value(const optional_t<T>& o)
{
    return o.impl.value();
}

template<typename T, typename F>
auto map(const optional_t<T>& o, F f)
{
    return has_value(o) ? optional(f(value(o))) : optional<std::invoke_result_t<F, T>>();
}

template<typename T, typename U>
auto zip(const optional_t<T>& o, const optional_t<U>& p)
{
    return has_value(o) && has_value(p) ? optional(std::pair(value(o), value(p))) : optional<std::pair<T, U>>();
}

template<typename T, typename U>
auto zip(const T& t, const optional_t<U>& p)
{
    return zip(optional(t), p);
}

template<typename T, typename U>
auto zip(const optional_t<T>& o, const U& u)
{
    return zip(o, optional(u));
}




//=============================================================================
template<typename T, typename U> auto operator+(T a, U b) { return map(zip(a, b), [] (auto t) { return std::apply(std::plus<>(), t); }); }
template<typename T, typename U> auto operator-(T a, U b) { return map(zip(a, b), [] (auto t) { return std::apply(std::minus<>(), t); }); }
template<typename T, typename U> auto operator*(T a, U b) { return map(zip(a, b), [] (auto t) { return std::apply(std::multiplies<>(), t); }); }
template<typename T, typename U> auto operator/(T a, U b) { return map(zip(a, b), [] (auto t) { return std::apply(std::divides<>(), t); }); }
template<typename T> auto operator+(optional_t<T> a) { return map(a, [] (auto x) { return x; }); }
template<typename T> auto operator-(optional_t<T> a) { return map(a, std::negate<>()); }

} // namespace numeric




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_numeric_optional()
{
    require(  has_value(numeric::optional(12)));
    require(! has_value(numeric::optional<int>()));
    require(! has_value(map(numeric::optional<double>(), std::negate<>())));
    require(  value(map(numeric::optional(1.0), std::negate<>())) == -1.0);
    require(  has_value(zip(numeric::optional(12), numeric::optional(13.0))));
    require(! has_value(zip(numeric::optional(12), numeric::optional<double>())));
    require(! has_value(zip(numeric::optional<int>(), numeric::optional(13.0))));
    require(! has_value(zip(numeric::optional<int>(), numeric::optional<double>())));
    require(  value(zip(12, numeric::optional(13.0))) == std::pair(12, 13.0));
    require(  value(zip(numeric::optional(12), 13.0)) == std::pair(12, 13.0));

    auto a = numeric::optional(12);
    auto b = numeric::optional(13.0);
    require(value(a + b) == 25.0);
    require(value(a + 13.0) == 25.0);
    require(value(+a) == +12);
    require(value(-a) == -12);
}

#endif // DO_UNIT_TESTS
