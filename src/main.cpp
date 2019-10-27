
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
#include "binary_serialize.hpp"
#include "bqo_tree.hpp"
#include "bsp_tree.hpp"
#include "dimensional.hpp"
#include "hdf5.hpp"
#include "linked_list.hpp"
#include "ndarray.hpp"
#include "numeric_array.hpp"
#include "numeric_optional.hpp"
#include "numeric_tuple.hpp"
#include "rational.hpp"
#include "sequence.hpp"
#include "unit_test.hpp"




//=============================================================================
int main()
{
    start_unit_tests();
    test_binary_serialize();
    test_bqo_tree();
    test_bsp_tree();
    test_dimensional();
    test_hdf5();
    test_linked_list();
    test_ndarray();
    test_numeric_array();
    test_numeric_optional();
    test_numeric_tuple();
    test_rational();
    test_sequence();
    report_test_results();
    return 0;
}
