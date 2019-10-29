
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




#define DO_UNIT_TESTS
#include "physics_euler.hpp"




//=============================================================================
int test_physics()
{
    auto require_cons_to_prim = [] (euler::primitive_t p0)
    {
        auto p1 = euler::recover_primitive(euler::conserved_density(p0, 1.4), 1.4);
        require(all(map(p1 - p0, [] (auto x) { return std::abs(x.value) < 1e-12; })));
    };

    require_cons_to_prim(euler::primitive(1.0, 1.0));
    require_cons_to_prim(euler::primitive(1.5, 1e3));
    require_cons_to_prim(euler::primitive(1.0, 0.0, 0.0, 0.0, 1.0));
    require_cons_to_prim(euler::primitive(5.7, 1.0, 2.0, 3.0, 5.0));

    return 0;
}
