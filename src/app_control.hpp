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
#include <chrono>
#include <tuple>
#include "core_sequence.hpp"




//=============================================================================
namespace control
{




using time_point_t = std::chrono::high_resolution_clock::time_point;




/**
 * A typedef intended for use with T as a solution state. The inner pair is a
 * solution state combined with a time point, helpful for use as aœ performance
 * diagnostic. The outer pair is two solution-states-time-point pairs.
 * Convenience functions below aid in accessing these four effective
 * "data members" by name.
 */
template<typename T>
using timed_pair_t = std::pair<std::pair<T, time_point_t>, std::pair<T, time_point_t>>;




template<typename T> auto last_time_point(timed_pair_t<T> p) { return p.first.second; }
template<typename T> auto this_time_point(timed_pair_t<T> p) { return p.second.second; }
template<typename T> auto last_state(timed_pair_t<T> p) { return p.first.first; }
template<typename T> auto this_state(timed_pair_t<T> p) { return p.second.first; }
template<typename T> auto microseconds_separating(timed_pair_t<T> p)
{
    auto t0 = last_time_point(p);
    auto t1 = this_time_point(p);
    return 1e-3 * std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
};
template<typename T> auto milliseconds_separating(timed_pair_t<T> p)
{
    return 1e3 * microseconds_separating(p);
};




/**
 * @brief      Return a sequence of high resolution time points
 *
 * @return     A sequence of std::high_resolution_clock::time_point
 */
inline auto time_point_sequence()
{
    using namespace std::chrono;
    return seq::generate(high_resolution_clock::now(), [] (auto) { return high_resolution_clock::now(); });
}

} // namespace control
