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
#include <cmath>
#include <algorithm>
#include "core_numeric_tuple.hpp"




//=============================================================================
namespace mara {




//=============================================================================
template<typename CoordinatesType>
struct plm_gradient_scalar_t
{
    template<typename T>
    static auto minabs(T a, T b, T c)
    {
        return std::min(std::min(abs(a), abs(b)), abs(c));
    }

    template<typename T>
    static auto sgn(T x)
    {
        return copysign(1.0, x);
    }

    template<typename T>
    auto operator()(numeric::tuple_t<T, T, T> y) const
    {
        auto [yl, y0, yr] = y.impl;
        auto a = (y0 - yl) / (x0 - xl) * theta;
        auto b = (yr - yl) / (xr - xl);
        auto c = (yr - y0) / (xr - x0) * theta;
        return 0.25 * abs(sgn(a) + sgn(b)) * (sgn(a) + sgn(c)) * minabs(a, b, c);
    }
    CoordinatesType xl, x0, xr;
    double theta = 1.5;
};




//=============================================================================
struct plm_gradient_t
{
    template<typename... Ts, std::size_t... Is>
    static auto zip_tuple_impl(numeric::tuple_t<Ts...> t, numeric::tuple_t<Ts...> u, numeric::tuple_t<Ts...> v, std::index_sequence<Is...>)
    {
        return numeric::tuple(numeric::tuple(numeric::get<Is>(t), numeric::get<Is>(u), numeric::get<Is>(v))...);
    }

    template<typename... Ts>
    static auto zip_tuples(numeric::tuple_t<Ts...> t, numeric::tuple_t<Ts...> u, numeric::tuple_t<Ts...> v)
    {
        return zip_tuple_impl(t, u, v, std::make_index_sequence<sizeof...(Ts)>());
    }

    template<typename T>
    auto operator()(std::tuple<T, T, T> states) const
    {
        auto [pl, p0, pr] = states;
        return map(zip_tuples(pl, p0, pr), plm_gradient_scalar_t<double>{0.0, 1.0, 2.0, theta});
    }

    template<typename X, typename P>
    auto operator()(std::tuple<std::tuple<X, X, X>, std::tuple<P, P, P>> coordinates_states) const
    {
        auto [coordinates, states] = coordinates_states;
        auto [xl, x0, xr] = coordinates;
        auto [pl, p0, pr] = states;
        return map(zip_tuples(pl, p0, pr), plm_gradient_scalar_t<X>{xl, x0, xr, theta});
    }
    double theta = 1.5;
};




//=============================================================================
inline auto plm_gradient(double theta)
{
    return plm_gradient_t{theta};
}

} // namespace mara
