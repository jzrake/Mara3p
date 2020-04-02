/**
 ==============================================================================
 Copyright 2020, Jonathan Zrake

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
#include <map>
#include <mutex>
#include <tuple>




/**
 * @brief      Drop-in replacement for std::invoke which memoizes the return
 *             result.
 *
 * @param[in]  function  The function whose results need to be cached
 * @param[in]  args      The function arguments
 *
 * @tparam     Function  The function type
 * @tparam     Args      The argument types
 *
 * @return     A value obtained either by evaluating the function, or by
 *             recalling it from a cache.
 *
 * @note       The function provided must not be a type-erase function object
 *             like a raw function pointer or std::function, because this
 *             function depends on the uniqueness of the Function template
 *             parameter. If you were to call invoke_memoized(f, a) and
 *             invoke_memoized(g, b) in the same translation unit, where f and g
 *             were function pointers of the same type, and a and b were
 *             arguments of the same type, you'd end up using the same cache for
 *             both functions f and g. A reasonable attempt is made to detect
 *             these misuse cases via static_assert.
 */
template<typename Function, typename... Args>
auto invoke_memoized(Function function, Args... args)
{
    using key_type   = std::tuple<Args...>;
    using value_type = std::invoke_result_t<Function, Args...>;

    static_assert(! std::is_same_v<Function, std::function<value_type(Args...)>>,
        "cannot memoize on std::function (use a lambda instead)");

    static_assert(! std::is_same_v<Function, value_type(*)(Args...)>,
        "cannot memoize on function pointer (use a lambda instead)");

    static std::mutex mutex;
    static std::map<key_type, value_type> cache;

    auto key  = std::tuple(args...);
    auto lock = std::lock_guard<std::mutex>(mutex);

    if (cache.count(key))
    {
        return cache[key];
    }
    return cache[key] = std::apply(function, key);
}




/**
 * @brief      Return a memoized version of a single-argument function object.
 *             The returned function is only safe if invoked by a single thread
 *             at a time.
 *
 * @param[in]  function  The function to memoize.
 *
 * @tparam     ArgType   The function argument type
 * @tparam     Function  The function type
 *
 * @return     A memoized function object
 */
template<typename ArgType, typename Function>
auto memoize_not_thread_safe(Function function)
{
    using result_type = std::invoke_result_t<Function, ArgType>;

    return [function, cache = std::make_shared<std::map<ArgType, result_type>>()] (const ArgType& arg)
    {
        if (cache->count(arg))
        {
            return cache->at(arg);
        }
        return cache->emplace(arg, std::invoke(function, arg)).first->second;
    };
}
