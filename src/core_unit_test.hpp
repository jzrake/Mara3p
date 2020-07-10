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
#include <cstdio>




//=============================================================================
#define White          "\033[0m"
#define Red            "\033[0;31m"
#define Green          "\033[0;32m"
#define Blue           "\033[0;34m"
#define Cyan           "\033[0;36m"
#define Yellow         "\033[0;33m"
#define Grey           "\033[1;30m"
#define LightGrey      "\033[0;37m"
#define BrightRed      "\033[1;31m"
#define BrightGreen    "\033[1;32m"
#define BrightWhite    "\033[1;37m"
#define BrightYellow   "\033[1;33m"




//=============================================================================
void start_unit_tests();
void report_test_results();
void increment_pass_count();
void increment_fail_count();




//=============================================================================
#define require(expr) \
do { \
    if ((expr)) { \
        std::printf("%sTest passed ... [%s:%d] %s%s\n", Green, __FILE__, __LINE__, #expr, White); \
        increment_pass_count(); \
    } \
    else { \
        std::printf("%sTest failed ... [%s:%d] %s%s\n", BrightRed, __FILE__, __LINE__, #expr, White); \
        increment_fail_count(); \
    } \
} while (false)

#define require_throws(expr) \
do { \
    try { \
        expr; \
        increment_fail_count(); \
        std::printf("%sTest failed ... [%s:%d] %s%s [did not throw] \n", BrightRed, __FILE__, __LINE__, #expr, White); \
    } \
    catch (...) { \
        increment_pass_count(); \
        std::printf("%sTest passed ... [%s:%d] %s%s [threw] \n", Green, __FILE__, __LINE__, #expr, White); \
    } \
} while (false)
