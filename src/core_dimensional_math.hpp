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
#include "core_dimensional.hpp"




//=============================================================================
namespace dimensional {




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
inline void test_dimensional_math()
{
    require(sqrt(dimensional::mass).powers == std::tuple(0.5, 0.0, 0.0));
    require(pow(dimensional::mass, rational::number<1, 2>()).powers == std::tuple(0.5, 0.0, 0.0));
}

#endif // DO_UNIT_TESTS
