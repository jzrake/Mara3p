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
#include <tuple>
#include <functional>
#include <memory>




//=============================================================================
namespace nd {

namespace detail {

template<typename... Ts, typename... Us, std::size_t... Is>
auto zip_tuple_impl(std::tuple<Ts...> t, std::tuple<Us...> u, std::index_sequence<Is...>)
{
    return std::tuple(std::pair(std::get<Is>(t), std::get<Is>(u))...);
}

template<typename... Ts, typename... Us>
auto zip_tuples(std::tuple<Ts...> t, std::tuple<Us...> u)
{
    return zip_tuple_impl(t, u, std::make_index_sequence<sizeof...(Ts)>());
}

template<typename FunctionType, typename... Ts>
auto map(std::tuple<Ts...> t, FunctionType f)
{
    return std::apply([f] (auto... ts) { return std::tuple(f(ts)...); }, t);
}

template<typename... Ts>
auto check_uniform(const char* message, Ts... ts)
{
    auto vs = {std::common_type_t<Ts...>(ts)...};
    auto v0 = *std::begin(vs);

    for (auto v : vs)
        if (v != v0)
            throw std::invalid_argument(message);

    return v0;
}

} // namespace detail




//=============================================================================
using uint = unsigned long;

template<uint Rank>
struct uivec_t
{
    uint& operator[](std::size_t n) { return __value[n]; }
    const uint& operator[](std::size_t n) const { return __value[n]; }
    bool operator==(uivec_t b) const { for (std::size_t i = 0; i < Rank; ++i) { if (__value[i] != b[i]) return false; } return true; }
    bool operator!=(uivec_t b) const { for (std::size_t i = 0; i < Rank; ++i) { if (__value[i] != b[i]) return true; } return false; }
    uint __value[Rank];
};

template<typename... Args>
auto uivec(Args... args)
{
    return uivec_t<sizeof...(Args)>{uint(args)...};
}

template<uint Rank>
constexpr std::size_t size(uivec_t<Rank>)
{
    return Rank;
}

template<uint Rank>
auto begin(const uivec_t<Rank>& t)
{
    return std::begin(t.__value);
}

template<uint Rank>
auto end(const uivec_t<Rank>& t)
{
    return std::end(t.__value);
}

template<uint Index, uint Rank>
const uint& get(const uivec_t<Rank>& vec)
{
    static_assert(Index < Rank, "nd::get (invalid index)");
    return vec[Index];
}

template<uint Rank>
auto to_tuple(uivec_t<Rank> t)
{
    return apply([] (auto... is) { return std::tuple(is...); }, t);
}

template<uint Rank>
uivec_t<Rank> replace(uivec_t<Rank> vec, uint axis, uint value)
{
    if (axis >= Rank)
        throw std::out_of_range("nd::replace (invalid axis)");

    vec[axis] = value;
    return vec;
}

template<uint Rank>
uivec_t<Rank - 1> remove(uivec_t<Rank> t, std::size_t axis)
{
    if (axis >= Rank)
        throw std::out_of_range("nd::remove (invalid axis)");

    auto u = uivec_t<Rank - 1>{};

    for (std::size_t i = 0, j = 0; i < size(t); ++i)
        if (i != axis)
            u[j++] = t[i];
    return u;
}

template<uint Rank>
uivec_t<Rank + 1> insert(uivec_t<Rank> t, std::size_t axis, uint value)
{
    if (axis > Rank)
        throw std::out_of_range("nd::insert (invalid axis)");

    auto u = uivec_t<Rank + 1>{};

    for (std::size_t i = 0, j = 0; i < size(u); ++i)
    {
        if (i != axis)
            u[i] = t[j++];
        else
            u[i] = value;
    }
    return u;
}

template<uint Rank>
uint product(uivec_t<Rank> t)
{
    uint n = 1;

    for (std::size_t i = 0; i < Rank; ++i)
        n *= t[i];
    return n;
}

template<uint Rank>
uint dot(uivec_t<Rank> t, uivec_t<Rank> u)
{
    uint n = 0;

    for (std::size_t i = 0; i < Rank; ++i)
        n += t[i] * u[i];
    return n;
}

template<uint Rank>
uivec_t<Rank> strides_row_major(uivec_t<Rank> shape)
{
    auto result = uivec_t<Rank>{};

    if constexpr (Rank > 0)
        result[Rank - 1] = 1;
    if constexpr (Rank > 1)
        for (int n = Rank - 2; n >= 0; --n)
            result[n] = result[n + 1] * shape[n + 1];

    return result;
}

template<uint Rank>
uivec_t<Rank> next(uivec_t<Rank> index, uivec_t<Rank> shape)
{
    if constexpr (Rank == 0)
        return index;

    auto n = Rank - 1;

    ++index[n];

    while (index[n] >= shape[n])
    {
        if (n == 0)
        {
            return shape;
        }
        index[n] = 0; --n; ++index[n];
    }
    return index;
}




//=============================================================================
namespace detail {

template<typename FunctionType, uint Rank, std::size_t... I>
constexpr decltype(auto) apply_impl(FunctionType&& f, uivec_t<Rank> t, std::index_sequence<I...>)
{
    return std::invoke(std::forward<FunctionType>(f), get<I>(t)...);
}

} // namespace detail

template<typename FunctionType, uint Rank>
constexpr decltype(auto) apply(FunctionType&& f, uivec_t<Rank> t)
{
    return detail::apply_impl(std::forward<FunctionType>(f), t, std::make_index_sequence<Rank>{});
}




//=============================================================================
template<typename ValueType>
class buffer_t
{
public:
    using value_type = ValueType;
    ~buffer_t() { delete [] memory; }
    buffer_t(std::size_t count) : memory(new value_type[count]) {}
    buffer_t(const buffer_t&) = delete;
    buffer_t& operator=(const buffer_t&) = delete;
    const value_type& operator[](std::size_t i) const { return memory[i]; }
    value_type& operator[](std::size_t i) { return memory[i]; }
    const value_type* data() const { return memory; }
    value_type* data() { return memory; }
private:
    value_type* memory;
};




//=============================================================================
template<typename ValueType, uint Rank>
struct shared_provider_t
{
    using value_type = ValueType;
    const value_type& operator()(uivec_t<Rank> i) const { return memory->operator[](dot(i, strides)); }
    const value_type* data() const { return memory ? memory->data() : nullptr; }

    std::shared_ptr<buffer_t<value_type>> memory;
    uivec_t<Rank> strides;
};




//=============================================================================
template<typename ValueType, uint Rank>
struct unique_provider_t
{
    using value_type = ValueType;
    const value_type& operator()(uivec_t<Rank> i) const { return memory->operator[](dot(i, strides)); }
    value_type& operator()(uivec_t<Rank> i)             { return memory->operator[](dot(i, strides)); }
    const value_type* data() const { return memory ? memory->data() : nullptr; }
    value_type* data()             { return memory ? memory->data() : nullptr; }
    auto shared() && { return shared_provider_t<value_type, Rank>{std::move(memory), strides}; }

    std::unique_ptr<buffer_t<value_type>> memory;
    uivec_t<Rank> strides;
};




//=============================================================================
template<typename ProviderType, uint Rank>
struct array_t
{
    using provider_type = ProviderType;
    using value_type = std::invoke_result_t<ProviderType, uivec_t<Rank>>;
    static constexpr uint rank = Rank;

    //=========================================================================
    struct iterator
    {
        using iterator_category = std::input_iterator_tag;
        using value_type        = uivec_t<Rank>;
        using difference_type   = std::ptrdiff_t;
        using pointer           = value_type*;
        using reference         = value_type&;

        iterator& operator++() { position = next(position, array.shape); return *this; }
        bool operator==(const iterator& other) const { return position == other.position; }
        bool operator!=(const iterator& other) const { return position != other.position; }
        decltype(auto) operator*() const { return array(position); }

        uivec_t<Rank> position;
        array_t array;
    };

    template<typename... Args>
    decltype(auto) operator()(Args... index_args) const { return provider(uivec(index_args...)); }
    decltype(auto) operator()(uivec_t<Rank> index) const { return provider(index); }
    decltype(auto) data() const { return provider.data(); }
    decltype(auto) data() { return provider.data(); }

    ProviderType provider;
    uivec_t<Rank> shape;
};




//=============================================================================
template<typename ValueType, uint Rank>
using shared_array = array_t<shared_provider_t<ValueType, Rank>, Rank>;

template<typename ValueType, uint Rank>
using unique_array = array_t<unique_provider_t<ValueType, Rank>, Rank>;

template<typename ValueType, uint Rank>
using dynamic_array = array_t<std::function<ValueType(uivec_t<Rank>)>, Rank>;




//=============================================================================
template<typename ProviderType, uint Rank>
auto make_array(ProviderType provider, uivec_t<Rank> shape)
{
    return array_t<ProviderType, Rank>{std::move(provider), shape};
}

template<typename ValueType, uint Rank>
auto make_unique_array(uivec_t<Rank> shape)
{
    auto t = strides_row_major(shape);
    auto p = nd::unique_provider_t<ValueType, Rank>{std::make_unique<nd::buffer_t<ValueType>>(product(shape)), t};
    return make_array(std::move(p), shape);    
}

template<typename ProviderType, uint Rank>
auto begin(const array_t<ProviderType, Rank>& array)
{
    return typename array_t<ProviderType, Rank>::iterator{{}, array};
}

template<typename ProviderType, uint Rank>
auto end(const array_t<ProviderType, Rank>& array)
{
    return typename array_t<ProviderType, Rank>::iterator{shape(array), array};
}

template<typename ProviderType, uint Rank>
constexpr uint rank(const array_t<ProviderType, Rank>&)
{
    return Rank;
}

template<typename ProviderType, uint Rank>
uint size(const array_t<ProviderType, Rank>& array)
{
    return product(array.shape);
}

template<typename ProviderType, uint Rank>
uivec_t<Rank> shape(const array_t<ProviderType, Rank>& array)
{
    return array.shape;
}

template<typename ProviderType, uint Rank>
uint shape(const array_t<ProviderType, Rank>& array, uint axis)
{
    if (axis >= Rank)
        throw std::out_of_range("nd::shape (invalid axis)");

    return array.shape[axis];
}




//=============================================================================
template<uint Rank>
struct index_space_row_major_t
{
    //=========================================================================
    struct iterator
    {
        using iterator_category = std::input_iterator_tag;
        using value_type        = uivec_t<Rank>;
        using difference_type   = std::ptrdiff_t;
        using pointer           = value_type*;
        using reference         = value_type&;

        iterator& operator++() { position = next(position, shape); return *this; }
        bool operator==(const iterator& other) const { return position == other.position; }
        bool operator!=(const iterator& other) const { return position != other.position; }
        auto operator*() const { return position; }

        uivec_t<Rank> position;
        uivec_t<Rank> shape;
    };
    uivec_t<Rank> shape;
};

template<uint Rank>
auto index_space(uivec_t<Rank> shape)
{
    return index_space_row_major_t<Rank>{shape};
}

template<uint Rank>
auto begin(const index_space_row_major_t<Rank>& space)
{
    return typename index_space_row_major_t<Rank>::iterator{{}, space.shape};
}

template<uint Rank>
auto end(const index_space_row_major_t<Rank>& space)
{
    return typename index_space_row_major_t<Rank>::iterator{space.shape, space.shape};
}




//=============================================================================
template<typename ValueType, uint Rank>
auto uniform(ValueType value, uivec_t<Rank> shape)
{
    return make_array([value] (const auto& i) { return value; }, shape);
}

template<typename ValueType=int, uint Rank>
auto zeros(uivec_t<Rank> shape)
{
    return uniform(ValueType(), shape);
}

template<typename ValueType=int, typename... Args>
auto zeros(Args... args)
{
    return uniform(ValueType(), uivec(args...));
}

template<typename... Args>
auto from(Args... args)
{
    std::common_type_t<Args...> a[] = {args...};

    return make_array(
        [a] (const auto& i) { return a[get<0>(i)]; },
        uivec(sizeof...(Args)));
}

inline auto range(long i0, long i1)
{
    if (i0 > i1)
        throw std::invalid_argument("nd::range (lower index larger than upper index)");

    return make_array([=] (const auto& i) { return i0 + get<0>(i); }, uivec(i1 - i0));
}

inline auto range(long i1)
{
    return range(0, i1);
}

inline auto linspace(double x0, double x1, uint count)
{
    return make_array(
        [=] (const auto& i) { return x0 + get<0>(i) * (x1 - x0) / (count - 1); },
        uivec(count));
}

template<typename... Args>
auto indexes(Args... args)
{
    return make_array([] (const auto& i) { return i; }, uivec(args...));
}

template<typename ProviderType, uint Rank>
auto indexes(const array_t<ProviderType, Rank>& array)
{
    return make_array([] (const auto& i) { return i; }, shape(array));
}




//=============================================================================
template<typename... ProviderTypes>
auto cartesian_product(array_t<ProviderTypes, 1>... arrays)
{
    auto f = [array_tuple = std::tuple(arrays...)] (const uivec_t<sizeof...(arrays)>& i)
    {
        return apply([] (auto&&... args)
        {
            return std::tuple(std::get<0>(args)(std::get<1>(args))...);
        }, detail::zip_tuples(array_tuple, to_tuple(i)));
    };
    return make_array(f, uivec(size(arrays)...));
}

template<typename... ProviderTypes, uint Rank>
auto zip(array_t<ProviderTypes, Rank>... arrays)
{
    auto s = detail::check_uniform("nd::zip (arrays have non-uniform shapes)", shape(arrays)...);
    auto f = [arrays...] (const auto& i)
    {
        return std::tuple(arrays(i)...);
    };
    return make_array(f, s);
}

template<std::size_t Index, typename ProviderType, uint Rank>
auto get(array_t<ProviderType, Rank> array)
{
    return map(array, [] (const auto& x) { return get<Index>(x); });
}

template<typename ProviderType, uint Rank, typename FunctionType>
auto map(array_t<ProviderType, Rank> array, FunctionType function)
{
    return make_array(
        [array, function] (const auto& i) { return function(array(i)); },
        shape(array));
}

template<typename ProviderType1, typename ProviderType2, uint Rank>
auto concat(array_t<ProviderType1, Rank> array1, array_t<ProviderType2, Rank> array2, uint axis)
{
    if (axis >= Rank)
        throw std::invalid_argument("nd::concat (axis must be smaller than array rank)");

    if (remove(shape(array1), axis) != remove(shape(array2), axis))
        throw std::invalid_argument("nd::concat (incompatible array shapes)");

    auto n1 = shape(array1, axis);
    auto n2 = shape(array2, axis);

    return make_array(
        [=] (const auto& i)
        {
            return i[axis] < n1 ? array1(i) : array2(replace(i, axis, i[axis] - n1));
        },
        replace(shape(array1), axis, n1 + n2));
}

template<typename ProviderType, uint Rank>
auto freeze(array_t<ProviderType, Rank> array, uint axis, uint index)
{
    if (axis >= Rank)
        throw std::invalid_argument("nd::freeze (axis is not smaller than array rank)");

    if (index >= shape(array)[axis])
        throw std::out_of_range("nd::freeze (index out of range)");

    return make_array(
        [=] (const auto& i) { return array(insert(i, axis, index)); },
        remove(shape(array), axis));
}

template<typename ProviderType, uint Rank>
auto new_axis(array_t<ProviderType, Rank> array, uint axis)
{
    if (axis > Rank)
        throw std::invalid_argument("nd::new_axis (axis is larger than array rank)");

    return make_array(
        [=] (const auto& i) { return array(remove(i, axis)); },
        insert(shape(array), axis, 1));
}

template<typename ProviderType, uint Rank>
auto select(array_t<ProviderType, Rank> array, uint axis, long start_s, long final_s)
{
    if (axis > Rank)
        throw std::invalid_argument("nd::select (axis is larger than array rank)");

    auto start = uint(start_s >= 0 ? start_s : shape(array, axis) + start_s);
    auto final = uint(final_s >= 0 ? final_s : shape(array, axis) + final_s);

    if (start > final)
        throw std::invalid_argument("nd::select (lower index larger than upper index)");

    if (start >= shape(array, axis) || final > shape(array, axis))
        throw std::out_of_range("nd::select (selection out of range)");

    return make_array(
        [=] (const auto& i) { return array(replace(i, axis, i[axis] + start)); },
        replace(shape(array), axis, final - start));
}

template<typename ProviderType, uint Rank>
auto select(array_t<ProviderType, Rank> array, uint axis, long start_s)
{
    if (axis > Rank)
        throw std::invalid_argument("nd::select (axis is larger than array rank)");

    return select(array, axis, start_s, shape(array, axis));
}

template<typename ProviderType, uint Rank>
auto to_shared(array_t<ProviderType, Rank> array)
{
    using value_type = typename array_t<ProviderType, Rank>::value_type;

    auto s = shape(array);
    auto t = strides_row_major(s);
    auto p = unique_provider_t<value_type, Rank>{std::make_unique<buffer_t<value_type>>(size(array)), t};

    for (auto i : index_space(s))
        p(i) = array(i);

    return make_array(std::move(p).shared(), s);
}

template<typename ProviderType, uint Rank>
auto to_dynamic(array_t<ProviderType, Rank> array)
{
    using value_type = typename array_t<ProviderType, Rank>::value_type;
    return make_array(std::function<value_type(uivec_t<Rank>)>(array.provider), shape(array));
}




//=============================================================================
template<std::size_t Index>
auto get()
{
    return [] (auto a) { return get<Index>(a); };
}

template<typename FunctionType>
auto map(FunctionType f)
{
    return [f] (auto a) { return map(a, f); };
}

template<typename ProviderType, uint Rank>
inline auto concat(array_t<ProviderType, Rank> b, uint axis)
{
    return [=] (auto a) { return concat(a, b, axis); };
}

inline auto freeze    (uint a, uint i)         { return [=] (auto array) { return freeze     (array, a, i); }; }
inline auto new_axis  (uint a)                 { return [=] (auto array) { return new_axis   (array, a); }; }
inline auto select    (uint a, long s, long f) { return [=] (auto array) { return select     (array, s, f); }; }
inline auto select    (uint a, long s)         { return [=] (auto array) { return select     (array, s); }; }
inline auto to_shared ()                       { return [ ] (auto array) { return to_shared  (array); }; }
inline auto to_dynamic()                       { return [ ] (auto array) { return to_dynamic (array); }; }




//=============================================================================
template<typename ProviderType, uint Rank>
auto where(array_t<ProviderType, Rank> array)
{
    auto b = map(array, [] (auto x) { return int(bool(x)); });
    auto p = unique_provider_t<uivec_t<Rank>, Rank>{std::make_unique<buffer_t<uivec_t<Rank>>>(sum(b)), {1}};
    auto n = uint(0);

    for (auto i : indexes(array))
        if (array(i))
            p(uivec(n++)) = i;

    return make_array(std::move(p).shared(), uivec(n));
}

template<typename ProviderType, uint Rank, typename FunctionType>
auto reduce(array_t<ProviderType, Rank> array, FunctionType function, typename array_t<ProviderType, Rank>::value_type x)
{
    for (auto y : array)
        x = function(x, y);
    return x;
}

template<typename ProviderType, uint Rank>
auto sum(array_t<ProviderType, Rank> array)
{
    return reduce(array, std::plus<>(), 0);
}

template<typename ProviderType, uint Rank>
auto product(array_t<ProviderType, Rank> array)
{
    return reduce(array, std::multiplies<>(), 1);
}

template<typename ProviderType, uint Rank>
bool any(array_t<ProviderType, Rank> array)
{
    for (auto x : array)
        if (x)
            return true;
    return false;
}

template<typename ProviderType, uint Rank>
bool all(array_t<ProviderType, Rank> array)
{
    for (auto x : array)
        if (! x)
            return false;
    return true;
}

inline auto sum()     { return [] (auto array) { return sum    (array); }; }
inline auto product() { return [] (auto array) { return product(array); }; }
inline auto any()     { return [] (auto array) { return any    (array); }; }
inline auto all()     { return [] (auto array) { return all    (array); }; }




//=============================================================================
template<typename A, typename B, uint R> auto operator+ (array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::plus<>(),            ab); }); }
template<typename A, typename B, uint R> auto operator- (array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::minus<>(),           ab); }); }
template<typename A, typename B, uint R> auto operator* (array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::multiplies<>(),      ab); }); }
template<typename A, typename B, uint R> auto operator/ (array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::divides<>(),         ab); }); }
template<typename A, typename B, uint R> auto operator< (array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::less<>(),            ab); }); }
template<typename A, typename B, uint R> auto operator> (array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::greater<>(),         ab); }); }
template<typename A, typename B, uint R> auto operator<=(array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::less_equal<>(),      ab); }); }
template<typename A, typename B, uint R> auto operator>=(array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::greater_equal<>(),   ab); }); }
template<typename A, typename B, uint R> auto operator==(array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::equal_to<>(),        ab); }); }
template<typename A, typename B, uint R> auto operator!=(array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::not_equal_to<>(),    ab); }); }
template<typename A, typename B, uint R> auto operator&&(array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::logical_and<>(),     ab); }); }
template<typename A, typename B, uint R> auto operator||(array_t<A, R> a, array_t<B, R> b) { return map(zip(a, b), [] (auto ab) { return apply(std::logical_or<>(),      ab); }); }
template<typename A, typename T, uint R> auto operator+ (array_t<A, R> a, T b) { return a +  uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator- (array_t<A, R> a, T b) { return a -  uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator* (array_t<A, R> a, T b) { return a *  uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator/ (array_t<A, R> a, T b) { return a /  uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator< (array_t<A, R> a, T b) { return a <  uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator> (array_t<A, R> a, T b) { return a >  uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator<=(array_t<A, R> a, T b) { return a <= uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator>=(array_t<A, R> a, T b) { return a >= uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator==(array_t<A, R> a, T b) { return a == uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator!=(array_t<A, R> a, T b) { return a != uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator&&(array_t<A, R> a, T b) { return a && uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator||(array_t<A, R> a, T b) { return a || uniform(b, shape(a)); }
template<typename A, typename T, uint R> auto operator+ (T a, array_t<A, R> b) { return uniform(a, shape(b)) +  b; }
template<typename A, typename T, uint R> auto operator- (T a, array_t<A, R> b) { return uniform(a, shape(b)) -  b; }
template<typename A, typename T, uint R> auto operator* (T a, array_t<A, R> b) { return uniform(a, shape(b)) *  b; }
template<typename A, typename T, uint R> auto operator/ (T a, array_t<A, R> b) { return uniform(a, shape(b)) /  b; }
template<typename A, typename T, uint R> auto operator< (T a, array_t<A, R> b) { return uniform(a, shape(b)) >  b; }
template<typename A, typename T, uint R> auto operator> (T a, array_t<A, R> b) { return uniform(a, shape(b)) <  b; }
template<typename A, typename T, uint R> auto operator<=(T a, array_t<A, R> b) { return uniform(a, shape(b)) <= b; }
template<typename A, typename T, uint R> auto operator>=(T a, array_t<A, R> b) { return uniform(a, shape(b)) >= b; }
template<typename A, typename T, uint R> auto operator==(T a, array_t<A, R> b) { return uniform(a, shape(b)) == b; }
template<typename A, typename T, uint R> auto operator!=(T a, array_t<A, R> b) { return uniform(a, shape(b)) != b; }
template<typename A, typename T, uint R> auto operator&&(T a, array_t<A, R> b) { return uniform(a, shape(b)) && b; }
template<typename A, typename T, uint R> auto operator||(T a, array_t<A, R> b) { return uniform(a, shape(b)) || b; }
template<typename A, uint R> auto operator+(array_t<A, R> a) { return a; };
template<typename A, uint R> auto operator-(array_t<A, R> a) { return map(a, std::negate<>()); };
template<typename A, uint R, typename F> auto operator|(array_t<A, R> a, F f) { return f(a); }

} // namespace nd




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_ndarray()
{
    // remove(uivec)
    require(remove(nd::uivec(0, 1, 2), 0) == nd::uivec(1, 2));
    require(remove(nd::uivec(0, 1, 2), 1) == nd::uivec(0, 2));
    require(remove(nd::uivec(0, 1, 2), 2) == nd::uivec(0, 1));

    // insert(uivec)
    require(insert(nd::uivec(0, 1, 2), 0, 9) == nd::uivec(9, 0, 1, 2));
    require(insert(nd::uivec(0, 1, 2), 1, 9) == nd::uivec(0, 9, 1, 2));
    require(insert(nd::uivec(0, 1, 2), 2, 9) == nd::uivec(0, 1, 9, 2));
    require(insert(nd::uivec(0, 1, 2), 3, 9) == nd::uivec(0, 1, 2, 9));

    // factories: zeros, range, from, linspace
    require(all(nd::zeros<double>(3, 4, 5) == nd::zeros(nd::uivec(3, 4, 5))));
    require(all(nd::range(4, 7) == nd::from(4, 5, 6)));
    require(all(nd::linspace(0, 1, 10) == nd::linspace(0, 2, 10) * 0.5));
    require_throws(nd::range(1, 0));

    // all, any, operator==, operator!=
    require(all(nd::linspace(0, 1, 10) == nd::linspace(0, 1, 10)) == true);
    require(any(nd::linspace(0, 1, 10) != nd::linspace(0, 1, 10)) == false);

    // to_shared, to_dynamic
    static_assert(std::is_same_v<nd::shared_array <nd::uint, 1>, decltype(to_shared (nd::range(10)))>);
    static_assert(std::is_same_v<nd::dynamic_array<nd::uint, 1>, decltype(to_dynamic(nd::range(10)))>);
    require(all(nd::range(10) == to_shared(nd::range(10))));
    require(all(nd::range(10) == to_dynamic(nd::range(10))));
    require(*to_shared(nd::range(5, 10)).data() == 5);

    // indexes, freeze, get
    require(all(nd::get<0>(freeze(nd::indexes(3, 3), 1, 2)) == nd::from(0, 1, 2)));
    require(all(nd::get<1>(freeze(nd::indexes(3, 3), 1, 2)) == nd::from(2, 2, 2)));

    // concat
    static_assert(std::is_same_v<decltype(concat(nd::zeros<int>(3, 4), nd::zeros<double>(3, 4), 2))::value_type, double>);
    require(shape(concat(nd::zeros(3, 4, 5), nd::zeros(6, 4, 5), 0)) == nd::uivec(9, 4, 5));
    require(shape(concat(nd::zeros(3, 4, 5), nd::zeros(3, 5, 5), 1)) == nd::uivec(3, 9, 5));
    require(shape(concat(nd::zeros(3, 4, 5), nd::zeros(3, 4, 4), 2)) == nd::uivec(3, 4, 9));
    require(concat(nd::indexes(3, 4, 5), nd::indexes(3, 4, 5), 0)(0, 0, 0) == nd::indexes(3, 4, 5)(0, 0, 0));
    require(concat(nd::indexes(3, 4, 5), nd::indexes(3, 4, 5), 0)(3, 0, 0) == nd::indexes(3, 4, 5)(0, 0, 0));

    require_throws(concat(nd::zeros(3, 4, 5), nd::zeros(3, 4, 5), 3)); // invalid axis
    require_throws(concat(nd::zeros(3, 4, 5), nd::zeros(4, 4, 5), 2)); // mismatch on axes 0, 1
    require_throws(concat(nd::zeros(3, 4, 5), nd::zeros(3, 4, 4), 0)); // mismatch on axes 1, 2
    require_throws(concat(nd::zeros(3, 4, 5), nd::zeros(3, 4, 4), 1)); // mismatch on axes 2, 0

    // select
    require(size(select(nd::linspace(0, 1, 10), 0, 5, 10)) == 5);
    require(select(nd::linspace(0, 1, 10), 0, 5, 10)(0) == nd::linspace(0, 1, 10)(5));
    require(select(nd::linspace(0, 1, 10), 0, 5, 10)(4) == nd::linspace(0, 1, 10)(9));
    require(all(select(nd::indexes(2, 3), 0, 0, 2) == nd::indexes(2, 3)));
    require(all(select(nd::indexes(2, 3), 0, 0, 1) == nd::indexes(1, 3)));
    require(all(select(nd::range(10), 0, -1) == nd::from(9)));
    require(all(select(nd::range(10), 0, 0, -9) == nd::from(0)));
    require(all(select(nd::range(10), 0, 3, -3) == nd::from(3, 4, 5, 6)));

    require_throws(select(nd::range(10), 0, 2, 1));
    require_throws(select(nd::range(10), 0, 0, 11));
    require_throws(select(nd::range(10), 1, 0, 1));

    // sum, product
    require(sum(nd::range(1, 4)) == 6);
    require(product(nd::range(1, 4)) == 6);

    // where, operator<, operator>, operator<=, operator>=
    require(size(where(nd::range(10) < 5)) == 5);
    require(size(where(nd::range(10) > 5)) == 4);
    require(size(where(nd::range(10) <= 5)) == 6);
    require(size(where(nd::range(10) >= 5)) == 5);

    // operator|
    require(all((nd::range(2) | nd::map([] (auto i) { return 2 * i; })) == nd::from(0, 2)));

    // operator+, operator-
    require(all(+nd::range(10) == nd::range(10)));
    require(all(-nd::range(10) == 0 - nd::range(10)));
}

#endif // DO_UNIT_TESTS
