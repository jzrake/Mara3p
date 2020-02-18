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
#include <functional>
#include <string>
#include <tuple>
#include "core_dimensional.hpp"
#include "core_rational.hpp"




//=============================================================================
namespace control
{




//=============================================================================
struct task_t
{
    std::string name;
    dimensional::unit_time next_time = 0.0;
    unsigned long count = 0;
};




//=============================================================================
inline task_t task(std::string name, dimensional::unit_time next_time=0.0)
{
    return {name, next_time, 0};
}

inline auto jump(task_t task, dimensional::unit_time time, dimensional::unit_time cadence)
{
    if (time >= task.next_time)
    {
        task.count += 1;
        task.next_time = task.next_time + cadence;
    }
    return task;
}

inline auto jump(task_t task, dimensional::unit_time time, double factor, dimensional::unit_time reference_time=0.0)
{
    if (time >= task.next_time)
    {
        task.count += 1;
        task.next_time = (task.next_time - reference_time) * factor + reference_time;
    }
    return task;
}

inline auto jump_task(dimensional::unit_time time, std::function<dimensional::unit_time(std::string)> cadence)
{
    return [=] (task_t task)
    {
        return jump(task, time, cadence(task.name));
    };
}





using time_point_t = std::chrono::high_resolution_clock::time_point;





/**
 * A typedef intended for use with T as a solution state. The inner pair is a
 * solution state combined with a time point, helpful for use as a performance
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

template<typename Function, typename... Args>
auto invoke_timed(Function f, Args... args)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    auto result = std::invoke(f, args...);
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::pair(result, t1 - t0);
}




/**
 * @brief      Advance a solution state using a common low-storage explicit
 *             Runge-Kutta scheme of order 1, 2, or 3.
 *
 * @param[in]  u_plus_du          The base scheme (single step)
 * @param[in]  rk_order           The order
 * @param[in]  u0                 The solution state to step forward from
 *
 * @tparam     BaseSchemeType     The base scheme type
 * @tparam     SolutionStateType  The solutions state type
 *
 * @return     An updated solution state
 *
 * @note       SolutionStateType must define an overload of the weighted_sum
 *             function. For two states s0 and s1 this function must return b *
 *             s0 + (1 - b) * s1, for a rational number b.
 */
template<typename BaseSchemeType, typename SolutionStateType>
auto advance_runge_kutta(BaseSchemeType u_plus_du, unsigned rk_order, SolutionStateType u0)
{
    switch (rk_order)
    {
        case 1:
        {
            return u_plus_du(u0);
        }
        case 2:
        {
            auto u1 = u0;
            u1 = u_plus_du(u0);
            u1 = weighted_sum(u0, u_plus_du(u1), rational::number(1, 2));
            return u1;
        }
        case 3:
        {
            auto u1 = u0;
            u1 = u_plus_du(u0);
            u1 = weighted_sum(u0, u_plus_du(u1), rational::number(3, 4));
            u1 = weighted_sum(u0, u_plus_du(u1), rational::number(1, 3));
            return u1;
        }
    }
    throw std::invalid_argument("sedov::advance (rk_order must be 1, 2, or 3)");
}

template<typename BaseSchemeType>
auto advance_runge_kutta(BaseSchemeType u_plus_du, unsigned rk_order)
{
    return [u_plus_du, rk_order] (auto u0)
    {
        return advance_runge_kutta(u_plus_du, rk_order, u0);
    };
}

} // namespace control
