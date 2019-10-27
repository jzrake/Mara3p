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
#include "core_rational.hpp"




//=============================================================================
namespace dimensional {




template<long N, unsigned long D>
using static_rational_t = rational::static_number_t<N, D>;




//=============================================================================
template<long N1, long N2, long N3, unsigned long D1=1, unsigned long D2=1, unsigned long D3=1>
struct dimensions_t
{
    static_rational_t<N1, D1> mass;
    static_rational_t<N2, D2> dist;
    static_rational_t<N3, D3> time;
};

template<long N1, long N2, long N3, unsigned long D1=1, unsigned long D2=1, unsigned long D3=1>
struct quantity_t
{
    operator double()
    {
        static_assert(N1 == 0 && N2 == 0 && N3 == 0, "only scalars are implicitly convertable to double");
        return value;
    }
    double value = 0.0;
    static constexpr auto dimensions = dimensions_t<N1, N2, N3, D1, D2, D3>{};
    static constexpr auto powers = std::make_tuple(double(N1) / D1, double(N2) / D2, double(N3) / D3);
};




//=============================================================================
static quantity_t<1, 0, 0> mass {1.0};
static quantity_t<0, 1, 0> dist {1.0};
static quantity_t<0, 0, 1> time {1.0};




//=============================================================================
template<long N1, long N2, long N3, unsigned long D1=1, unsigned long D2=1, unsigned long D3=1>
auto dimensions(static_rational_t<N1, D1>, static_rational_t<N2, D2>, static_rational_t<N3, D3>)
{
    return dimensions_t<N1, N2, N3, D1, D2, D3>{};
}

template<long N1, long N2, long N3, unsigned long D1=1, unsigned long D2=1, unsigned long D3=1>
auto quantity(double value, dimensions_t<N1, N2, N3, D1, D2, D3>)
{
    return quantity_t<N1, N2, N3, D1, D2, D3>{value};
}

inline auto scalar(double value)
{
    return quantity_t<0, 0, 0>{value};
}




//=============================================================================
template<
long N1a, long N2a, long N3a, unsigned long D1a, unsigned long D2a, unsigned long D3a,
long N1b, long N2b, long N3b, unsigned long D1b, unsigned long D2b, unsigned long D3b>
auto operator*(dimensions_t<N1a, N2a, N3a, D1a, D2a, D3a> a, dimensions_t<N1b, N2b, N3b, D1b, D2b, D3b> b)
{
    return dimensions(a.mass + b.mass, a.dist + b.dist, a.time + b.time);
}

template<
long N1a, long N2a, long N3a, unsigned long D1a, unsigned long D2a, unsigned long D3a,
long N1b, long N2b, long N3b, unsigned long D1b, unsigned long D2b, unsigned long D3b>
auto operator/(dimensions_t<N1a, N2a, N3a, D1a, D2a, D3a> a, dimensions_t<N1b, N2b, N3b, D1b, D2b, D3b> b)
{
    return dimensions(a.mass - b.mass, a.dist - b.dist, a.time - b.time);
}

template<
long N, unsigned long D=1, long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto pow(dimensions_t<N1, N2, N3, D1, D2, D3> a, static_rational_t<N, D> p)
{
    return dimensions(a.mass * p, a.dist * p, a.time * p);
}

template<
long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto sqrt(dimensions_t<N1, N2, N3, D1, D2, D3> a)
{
    return pow(a, static_rational_t<1, 2>{});
}




//=============================================================================
template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto operator+(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return quantity_t<N1, N2, N3, D1, D2, D3>{a.value + b.value};
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto operator-(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return quantity_t<N1, N2, N3, D1, D2, D3>{a.value - b.value};
}

template<
long N1a, long N2a, long N3a, unsigned long D1a, unsigned long D2a, unsigned long D3a,
long N1b, long N2b, long N3b, unsigned long D1b, unsigned long D2b, unsigned long D3b>
auto operator*(quantity_t<N1a, N2a, N3a, D1a, D2a, D3a> a, quantity_t<N1b, N2b, N3b, D1b, D2b, D3b> b)
{
    return quantity(a.value * b.value, a.dimensions * b.dimensions);
}

template<
long N1a, long N2a, long N3a, unsigned long D1a, unsigned long D2a, unsigned long D3a,
long N1b, long N2b, long N3b, unsigned long D1b, unsigned long D2b, unsigned long D3b>
auto operator/(quantity_t<N1a, N2a, N3a, D1a, D2a, D3a> a, quantity_t<N1b, N2b, N3b, D1b, D2b, D3b> b)
{
    return quantity(a.value / b.value, a.dimensions / b.dimensions);
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto operator*(quantity_t<N1, N2, N3, D1, D2, D3> a, double s)
{
    return a * scalar(s);
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto operator*(double s, quantity_t<N1, N2, N3, D1, D2, D3> a)
{
    return scalar(s) * a;
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto operator/(quantity_t<N1, N2, N3, D1, D2, D3> a, double s)
{
    return a / scalar(s);
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto operator/(double s, quantity_t<N1, N2, N3, D1, D2, D3> a)
{
    return scalar(s) / a;
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
bool operator==(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return a.value == b.value;
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
bool operator!=(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return a.value != b.value;
}




//=============================================================================
template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto sqrt(quantity_t<N1, N2, N3, D1, D2, D3> a)
{
    return quantity(std::sqrt(a.value), pow(a.dimensions, static_rational_t<1, 2>{}));
}

template<long N, unsigned long D=1, long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
auto pow(quantity_t<N1, N2, N3, D1, D2, D3> a, static_rational_t<N, D> p={})
{
    return quantity(std::pow(a.value, p.value), pow(a.dimensions, p));
}

} // namespace dimensional




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
void test_dimensional()
{
    require(sizeof(dimensional::time) == sizeof(double));
    require((dimensional::mass * dimensional::mass).powers == std::make_tuple(2., 0., 0.));
    require((dimensional::mass / dimensional::mass).powers == std::make_tuple(0., 0., 0.));
    require((dimensional::dist * dimensional::dist).powers == std::make_tuple(0., 2., 0.));
    require((dimensional::dist / dimensional::dist).powers == std::make_tuple(0., 0., 0.));
    require((dimensional::time * dimensional::time).powers == std::make_tuple(0., 0., 2.));
    require((dimensional::time / dimensional::time).powers == std::make_tuple(0., 0., 0.));
    require((1. / dimensional::mass).powers == std::make_tuple(-1.,  0.,  0.));
    require((1. / dimensional::dist).powers == std::make_tuple( 0., -1.,  0.));
    require((1. / dimensional::time).powers == std::make_tuple( 0.,  0., -1.));
    require((2. * dimensional::mass).powers == std::make_tuple( 1.,  0.,  0.));
    require((2. * dimensional::dist).powers == std::make_tuple( 0.,  1.,  0.));
    require((2. * dimensional::time).powers == std::make_tuple( 0.,  0.,  1.));
    require(sqrt(dimensional::mass).powers == std::make_tuple(0.5, 0.0, 0.0));
    require(pow(dimensional::mass, rational::number<1, 2>()).powers == std::make_tuple(0.5, 0.0, 0.0));
}

#endif // DO_UNIT_TESTS
