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

auto extend_zeros(uint axis=0, uint count=1)
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

} // namespace nd
