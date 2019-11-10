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
#include "core_dimensional.hpp"
#include "core_ndarray.hpp"
#include "core_numeric_tuple.hpp"
#include "core_rational.hpp"




//=============================================================================
namespace mara {




//=============================================================================
template<typename ConservedType>
struct state_with_vertices_t
{
    rational::number_t iteration = 0;
    dimensional::unit_time time = 0.0;
    nd::shared_array<dimensional::unit_length, 1> vertices;
    nd::shared_array<ConservedType, 1> conserved;
};

template<typename ConservedType>
auto weighted_sum(state_with_vertices_t<ConservedType> s, state_with_vertices_t<ConservedType> t, rational::number_t b)
{
    return state_with_vertices_t<ConservedType>{
        s.iteration *         b + t.iteration *       (1 - b),
        s.time      * double(b) + t.time      * double(1 - b),
        s.vertices  * double(b) + t.vertices  * double(1 - b) | nd::to_shared(),
        s.conserved * double(b) + t.conserved * double(1 - b) | nd::to_shared(),
    };
}

} // namespace mara
