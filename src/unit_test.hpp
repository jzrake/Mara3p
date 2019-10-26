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




#pragma once
#include <iostream>




//=============================================================================
#define White          "[0m"
#define Red            "[0;31m"
#define Green          "[0;32m"
#define Blue           "[0;34m"
#define Cyan           "[0;36m"
#define Yellow         "[0;33m"
#define Grey           "[1;30m"
#define LightGrey      "[0;37m"
#define BrightRed      "[1;31m"
#define BrightGreen    "[1;32m"
#define BrightWhite    "[1;37m"
#define BrightYellow   "[1;33m"




static unsigned NumPassed = 0;
static unsigned NumFailed = 0;




//=============================================================================
#define require(expr) \
do { \
    if (expr) { \
        std::cout << '\033' << Green     << "Test passed ... [" << __FILE__ << "] " << #expr << std::endl; \
        ++NumPassed; \
    } \
    else { \
        std::cout << '\033' << BrightRed << "Test failed ... " << #expr << " [" << __FILE__ << " on line " << __LINE__ << "]" << std::endl; \
        ++NumFailed; \
    } \
} while (false)

#define require_throws(expr) \
do { \
    try { \
        expr; \
        ++NumFailed; \
        std::cout << '\033' << BrightRed << "Test failed ... " << #expr << " " << __FILE__ << " on line " << __LINE__ << " [did not throw]" << std::endl; \
    } \
    catch (...) { \
        ++NumPassed; \
        std::cout << '\033' << Green     << "Test passed ... " << #expr << " [threw]" << std::endl; \
    } \
} while (false)




//=============================================================================
inline void start_unit_tests()
{
    NumPassed = 0;
    NumFailed = 0;
}

inline void report_test_results()
{
    std::cout << std::endl << "===============================================================================" << std::endl;
    std::cout << NumPassed << " passed" << std::endl;

    if (NumFailed)
        std::cout << '\033' << BrightRed << NumFailed << " failed" << std::endl;        
}
