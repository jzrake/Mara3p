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
 LIABILITY, WHETHER IN AN ACTION OF CONT (eRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 ==============================================================================
*/




#include "core_unit_test.hpp"




static unsigned NumPassed = 0;
static unsigned NumFailed = 0;




//=============================================================================
void start_unit_tests()
{
    NumPassed = 0;
    NumFailed = 0;
}

void report_test_results()
{
    std::cout << std::endl << "===============================================================================" << std::endl;
    std::cout << NumPassed << " passed" << std::endl;

    if (NumFailed)
        std::cout << '\033' << BrightRed << NumFailed << " failed" << std::endl;        
}

void increment_pass_count()
{
    ++NumPassed;
}

void increment_fail_count()
{
    ++NumFailed;
}
