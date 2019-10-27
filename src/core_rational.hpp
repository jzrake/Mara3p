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
#include <numeric>




//=============================================================================
namespace rational {




//=============================================================================
struct number_t
{
    number_t(int value) : num(value), den(1) {}
    number_t(long num, unsigned long den) : num(num), den(den) {}
    operator double() { return double(num) / den; }

    long num = 0;
    unsigned long den = 1;
};

inline number_t number(long num, unsigned long den=1) { return {num, den}; }




//=============================================================================
inline number_t reduce(number_t a)
{
  return {a.num / long(std::gcd(std::abs(a.num), a.den)), a.den / std::gcd(std::abs(a.num), a.den)};
}

inline number_t reciprocal(number_t a)
{
    auto num = long(a.den) * ((a.num > 0) - (a.num < 0));
    auto den = (unsigned long)(std::abs(a.num));
    return {num, den};
}




//=============================================================================
inline number_t operator+(number_t a, number_t b) { return reduce({a.num * long(b.den) + long(a.den) * b.num, a.den * b.den}); }
inline number_t operator-(number_t a, number_t b) { return reduce({a.num * long(b.den) - long(a.den) * b.num, a.den * b.den}); }
inline number_t operator*(number_t a, number_t b) { return reduce({a.num * b.num, a.den * b.den}); }
inline number_t operator/(number_t a, number_t b) { return a * reciprocal(b); }

inline number_t& operator+=(number_t& a, number_t b) { a = a + b; return a; }
inline number_t& operator-=(number_t& a, number_t b) { a = a - b; return a; }
inline number_t& operator*=(number_t& a, number_t b) { a = a * b; return a; }
inline number_t& operator/=(number_t& a, number_t b) { a = a / b; return a; }

inline bool operator==(number_t a, number_t b) { return reduce(a / b).num == 1 && reduce(a / b).den == 1; }
inline bool operator!=(number_t a, number_t b) { return reduce(a / b).num != 1 || reduce(a / b).den != 1; }




//=============================================================================
template<long N, unsigned long D>
struct static_number_t
{
    static constexpr double value = double(N) / D;
    operator number_t() const { return number_t(N, D); }
};




//=============================================================================
template<long N, unsigned long D>
auto number()
{
    return static_number_t<N, D>{};
}

template<long N, unsigned long D>
auto reduce(static_number_t<N, D>)
{
    return static_number_t<N / long(std::gcd(N * ((N > 0) - (N < 0)), D)), D / std::gcd(N * ((N > 0) - (N < 0)), D)>{};
}

template<long N, unsigned long D>
auto reciprocal(static_number_t<N, D>)
{
    return static_number_t<long(D) * ((N > 0) - (N < 0)), (unsigned long)(N >= 0 ? N : -N)>{};
}

template<long N1, unsigned long D1, long N2, unsigned long D2>
auto operator+(static_number_t<N1, D1> a, static_number_t<N2, D2> b)
{
    return reduce(static_number_t<N1 * long(D2) + N2 * long(D1), D1 * D2>{});
}

template<long N1, unsigned long D1, long N2, unsigned long D2>
auto operator-(static_number_t<N1, D1> a, static_number_t<N2, D2> b)
{
    return reduce(static_number_t<N1 * long(D2) - N2 * long(D1), D1 * D2>{});
}

template<long N1, unsigned long D1, long N2, unsigned long D2>
auto operator*(static_number_t<N1, D1> a, static_number_t<N2, D2> b)
{
    return reduce(static_number_t<N1 * N2, D1 * D2>{});
}

template<long N1, unsigned long D1, long N2, unsigned long D2>
auto operator/(static_number_t<N1, D1> a, static_number_t<N2, D2> b)
{
    return a * reciprocal(b);
}

} // namespace rational




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
void test_rational()
{
    require(rational::number(+1, 2) + rational::number(+3, 2) == rational::number(+2));
    require(rational::number(+1, 2) - rational::number(+3, 2) == rational::number(-1));
    require(rational::number(+1, 2) * rational::number(+3, 2) == rational::number(+3, 4));
    require(rational::number(-1, 2) * rational::number(+3, 2) == rational::number(-3, 4));
    require(rational::number(+1, 2) / rational::number(+3, 2) == rational::number(+1, 3));
    require(rational::number(-1, 2) / rational::number(+3, 2) == rational::number(-1, 3));
    require(rational::number(+1, 2) / rational::number(+3, 2) == rational::number(+1, 3));
    require(rational::number(+1, 2) / rational::number(-3, 2) == rational::number(-1, 3));
    require(reciprocal(rational::number(+3, 2)) == rational::number(+2, 3));
    require(reciprocal(rational::number(-3, 2)) == rational::number(-2, 3));
    require(reduce(rational::number(+1, 3)) == rational::number(+1, 3));
    require(reduce(rational::number(-1, 3)) == rational::number(-1, 3));
    require((rational::number<+1, 3>().value == +1. / 3));
    require((rational::number<-1, 3>().value == -1. / 3));
    require((reduce(rational::number<+1, 3>()).value == +1. / 3));
    require((reduce(rational::number<-1, 3>()).value == -1. / 3));

    rational::static_number_t<3, 2> a, b;
    require((a + b).value == 3.0);
    require((a - b).value == 0.0);
    require((a * b).value == 9. / 4);
    require((a / b).value == 1.0);
}

#endif // DO_UNIT_TESTS
