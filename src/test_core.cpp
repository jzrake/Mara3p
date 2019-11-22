
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
#include "core_bqo_tree.hpp"
#include "core_bsp_tree.hpp"
#include "core_dimensional.hpp"
#include "core_dimensional_math.hpp"
#include "core_geometric.hpp"
#include "core_linked_list.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_numeric_array.hpp"
#include "core_numeric_optional.hpp"
#include "core_numeric_tuple.hpp"
#include "core_rational.hpp"
#include "core_sequence.hpp"
#include "core_unit_test.hpp"




//=============================================================================
int test_core()
{
    test_bqo_tree();
    test_bsp_tree();
    test_dimensional();
    test_dimensional_math();
    test_geometric();
    test_linked_list();
    test_ndarray();
    test_ndarray_ops();
    test_numeric_array();
    test_numeric_optional();
    test_numeric_tuple();
    test_rational();
    test_sequence();
    return 0;
}
