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
#include <string>




//=============================================================================
namespace util {




/**
 * @brief      Helper function to turn a function f(a, b, ...) of several
 *             variables into a function g(std::tuple(a, b, ...)) of a single
 *             tuple. This enables constructs like zip(A, B) | map(apply_to([]
 *             (auto a, auto b) { returh a + b; })), where zip turns functors A
 *             = F<T> and B = F<U> into a functor of tuples F<std::tuple<T, U>>.
 *
 * @param[in]  function      The function f(a, b, ...)
 *
 * @tparam     FunctionType  The type of the function f
 *
 * @return     A function g(std::tuple(a, b, ...)) of a single tuple
 */
template<typename FunctionType>
auto apply_to(FunctionType function)
{
    return [function] (auto t) { return std::apply(function, t); };
}




/**
 * @brief      Wrapper for the snprintf function.
 *
 * @param[in]  format_string  The c-style format string to use
 * @param[in]  args           The arguments for the format string
 *
 * @tparam     Args           The argument types
 *
 * @return     A formatted std::string.
 *
 * @note       This function is not meant to format long strings; the length of
 *             the result string is < 2048. Longer results are (safely)
 *             truncated.
 */
template<typename... Args>
std::string format(const char* format_string, Args... args)
{
    char result[2048];
    std::snprintf(result, 2048, format_string, args...);
    return result;
}

} // namespace util
