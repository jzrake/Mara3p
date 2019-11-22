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
#include "core_numeric_array.hpp"




//=============================================================================
namespace geometric
{




//=============================================================================
template<typename CoordinateType>
struct euclidean_vector_t
{
    CoordinateType component_1() const { return impl[0]; }
    CoordinateType component_2() const { return impl[1]; }
    CoordinateType component_3() const { return impl[2]; }
    CoordinateType component(unsigned axis) const { return impl.at(axis - 1); }

    numeric::array_t<CoordinateType, 3> impl;
};




//=============================================================================
struct unit_vector_t
{
    double component_1() const { return impl[0]; }
    double component_2() const { return impl[1]; }
    double component_3() const { return impl[2]; }
    double component(unsigned axis) const { return impl.at(axis - 1); }

    numeric::array_t<double, 3> impl;
};




//=============================================================================
template<typename T>
euclidean_vector_t<T> euclidean_vector(T x, T y, T z)
{
    return {numeric::array(x, y, z)};
}

template<typename T>
euclidean_vector_t<T> euclidean_vector(numeric::array_t<T, 3> v)
{
    return {v};
}

inline unit_vector_t unit_vector_on_axis(unsigned axis)
{
    switch (axis)
    {
        case 1: return unit_vector_t{{1.0, 0.0, 0.0}};
        case 2: return unit_vector_t{{0.0, 1.0, 0.0}};
        case 3: return unit_vector_t{{0.0, 0.0, 1.0}};
    }
    throw std::invalid_argument("geometric::unit_vector_on_axis (axis must be 1, 2, or 3)");
}




//=============================================================================
template<typename T, typename U>
auto dot(euclidean_vector_t<T> a, euclidean_vector_t<U> b)
{
    return sum(a.impl * b.impl);
}

template<typename T>
auto dot(unit_vector_t a, euclidean_vector_t<T> b)
{
    return sum(a.impl * b.impl);
}

template<typename T>
auto dot(euclidean_vector_t<T> b, unit_vector_t a)
{
    return sum(a.impl * b.impl);
}

template<typename T>
auto length_squared(euclidean_vector_t<T> a)
{
    return sum(a.impl * a.impl);
}

template<typename T, typename U>
auto cross(euclidean_vector_t<T> b, euclidean_vector_t<U> a)
{
    return euclidean_vector(a.impl[1] * b.impl[2], a.impl[2] * b.impl[0], a.impl[0] * b.impl[1]);
}




//=============================================================================
template<typename T, typename U> auto operator+(euclidean_vector_t<T> a, euclidean_vector_t<U> b) { return euclidean_vector(a.impl + b.impl); }
template<typename T, typename U> auto operator-(euclidean_vector_t<T> a, euclidean_vector_t<U> b) { return euclidean_vector(a.impl - b.impl); }
template<typename T, typename U> auto operator*(euclidean_vector_t<T> a, euclidean_vector_t<U> b) { return euclidean_vector(a.impl * b.impl); }
template<typename T, typename U> auto operator/(euclidean_vector_t<T> a, euclidean_vector_t<U> b) { return euclidean_vector(a.impl / b.impl); }
template<typename T, typename U> auto operator*(euclidean_vector_t<T> a, U b) { return euclidean_vector(a.impl * b); }
template<typename T, typename U> auto operator/(euclidean_vector_t<T> a, U b) { return euclidean_vector(a.impl / b); }
template<typename T, typename U> auto operator*(U b, euclidean_vector_t<T> a) { return euclidean_vector(b * a.impl); }
template<typename T, typename U> auto operator/(U b, euclidean_vector_t<T> a) { return euclidean_vector(b / a.impl); }
template<typename T, typename U> auto operator+(euclidean_vector_t<T> a) { return euclidean_vector(+a.impl); }
template<typename T, typename U> auto operator-(euclidean_vector_t<T> a) { return euclidean_vector(-a.impl); }
template<typename T, typename U> auto operator==(euclidean_vector_t<T> a, euclidean_vector_t<U> b) { return a.impl == b.impl; }
template<typename T, typename U> auto operator!=(euclidean_vector_t<T> a, euclidean_vector_t<U> b) { return a.impl != b.impl; }

template<typename T> auto operator*(unit_vector_t a, T b) { return euclidean_vector(a.impl * b); }
template<typename T> auto operator*(T a, unit_vector_t b) { return euclidean_vector(a * b.impl); }

} // namespace geometric




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_geometric()
{
    require(geometric::euclidean_vector(1, 2, 3).component_1() == 1);
    require(geometric::euclidean_vector(1, 2, 3).component_2() == 2);
    require(geometric::euclidean_vector(1, 2, 3) + geometric::euclidean_vector(1, 2, 3) == geometric::euclidean_vector(2, 4, 6));
    require(geometric::euclidean_vector(1, 2, 3) * 2.0 == geometric::euclidean_vector(2, 4, 6));
}

#endif // DO_UNIT_TESTS
