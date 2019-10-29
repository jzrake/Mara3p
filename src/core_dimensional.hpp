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
    quantity_t() {}
    quantity_t(double value) : value(value) {}

    operator double()
    {
        static_assert(N1 == 0 && N2 == 0 && N3 == 0, "only scalars are implicitly convertable to double");
        return value;
    }

    double value = 0.0;
    static constexpr auto dimensions = dimensions_t<N1, N2, N3, D1, D2, D3>{};
    static constexpr auto powers = std::tuple(double(N1) / D1, double(N2) / D2, double(N3) / D3);
};




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

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
bool operator<=(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return a.value <= b.value;
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
bool operator>=(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return a.value >= b.value;
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
bool operator<(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return a.value < b.value;
}

template<long N1, long N2, long N3, unsigned long D1, unsigned long D2, unsigned long D3>
bool operator>(quantity_t<N1, N2, N3, D1, D2, D3> a, quantity_t<N1, N2, N3, D1, D2, D3> b)
{
    return a.value > b.value;
}




//=============================================================================
using unit_mass             = quantity_t<1, 0, 0>;
using unit_length           = quantity_t<0, 1, 0>;
using unit_time             = quantity_t<0, 0, 1>;
using unit_velocity         = decltype(unit_length{}   / unit_time{});
using unit_acceleration     = decltype(unit_velocity{} / unit_time{});
using unit_momentum         = decltype(unit_velocity{} * unit_mass{});
using unit_force            = decltype(unit_momentum{} / unit_time{});
using unit_energy           = decltype(unit_force{}    * unit_length{});
using unit_area             = decltype(unit_length{}   * unit_length{});
using unit_volume           = decltype(unit_length{}   * unit_length{} * unit_length{});
using unit_mass_density     = decltype(unit_mass{}     / unit_volume{});
using unit_momentum_density = decltype(unit_momentum{} / unit_volume{});
using unit_energy_density   = decltype(unit_energy{}   / unit_volume{});
using unit_specific_energy  = decltype(unit_energy{}   / unit_mass{});

} // namespace dimensional




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_dimensional()
{
    auto mass = dimensional::unit_mass();
    auto dist = dimensional::unit_length();
    auto time = dimensional::unit_time();

    require(sizeof(time) == sizeof(double));
    require((mass * mass).powers == std::tuple(2., 0., 0.));
    require((mass / mass).powers == std::tuple(0., 0., 0.));
    require((dist * dist).powers == std::tuple(0., 2., 0.));
    require((dist / dist).powers == std::tuple(0., 0., 0.));
    require((time * time).powers == std::tuple(0., 0., 2.));
    require((time / time).powers == std::tuple(0., 0., 0.));
    require((1. / mass).powers == std::tuple(-1.,  0.,  0.));
    require((1. / dist).powers == std::tuple( 0., -1.,  0.));
    require((1. / time).powers == std::tuple( 0.,  0., -1.));
    require((2. * mass).powers == std::tuple( 1.,  0.,  0.));
    require((2. * dist).powers == std::tuple( 0.,  1.,  0.));
    require((2. * time).powers == std::tuple( 0.,  0.,  1.));
}

#endif // DO_UNIT_TESTS
