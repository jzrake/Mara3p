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
#include "app_config.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_rational.hpp"
#include "app_hdf5_std_map.hpp"
#include "app_hdf5_std_variant.hpp"
#include "app_serial.hpp"
#include "app_serial_ndarray.hpp"




//=============================================================================
int test_app()
{
    test_hdf5();
    test_hdf5_dimensional();
    test_hdf5_ndarray();
    test_hdf5_ndarray_dimensional();
    test_hdf5_numeric_array();
    test_hdf5_rational();
    test_hdf5_std_variant();
    test_hdf5_std_map();
    test_serial();
    test_serial_ndarray();
    return 0;
}
