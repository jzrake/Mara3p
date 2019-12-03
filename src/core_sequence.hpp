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
#include <any>
#include <array>
#include <memory>
#include <optional>




//=============================================================================
namespace seq {

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
auto map(const std::tuple<Ts...>& t, FunctionType fn)
{
    return std::apply([fn] (const auto&... ts) { return std::tuple(fn(ts)...); }, t);
}

template<typename... Ts>
auto values(const std::tuple<std::optional<Ts>...>& t)
{
    return std::apply([] (const auto&... os) { return std::tuple(os.value()...); }, t);
}

template<typename... Ts>
bool has_values(const std::tuple<std::optional<Ts>...>& t)
{
    auto impl = [] (const auto&... os) {
        return std::array{os.has_value()...};
    };
    for (auto has_value : std::apply(impl, t))
        if (! has_value)
            return false;
    return true;
}

template<typename... Ts>
auto optional_tuple(const std::tuple<std::optional<Ts>...>& t)
{
    return has_values(t) ? std::make_optional(values(t)) : std::optional<std::tuple<Ts...>>{};
}

template<typename T>
std::optional<std::any> optional_any(const std::optional<T>& o)
{
    return o.has_value() ? o.value() : std::optional<std::any>{};
}

template<typename T, typename U>
std::optional<std::pair<std::optional<T>, std::optional<U>>> optional_pair_of_either(const std::optional<T>& t, const std::optional<U>& u)
{
    return t.has_value() || u.has_value()
    ? std::optional<std::pair<std::optional<T>, std::optional<U>>>{{t, u}}
    : std::optional<std::pair<std::optional<T>, std::optional<U>>>{};
}

} // namespace detail




//=============================================================================
template<typename SequenceType>
struct position
{
    using type = typename decltype(start(std::declval<SequenceType>()))::value_type;
};

template<typename SequenceType>
using position_t = typename position<SequenceType>::type;

template<typename SequenceType>
struct value_type
{
    using type = decltype(obtain(std::declval<SequenceType>(), std::declval<position_t<SequenceType>>()));
};

template<typename SequenceType>
using value_type_t = typename value_type<SequenceType>::type;




//=============================================================================s
template<typename SequenceType>
struct iterator
{
    using sequence_type = SequenceType;
    using position_type = position_t<sequence_type>;
    using iterator_category = std::input_iterator_tag;
    using value_type        = value_type_t<SequenceType>;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type*;
    using reference         = value_type&;

    iterator& operator++() { position = next(sequence, position.value()); return *this; }
    bool operator!=(const iterator& other) const { return position.has_value() != other.position.has_value(); }
    auto operator*() const { return obtain(sequence, position.value()); }
    sequence_type sequence;
    std::optional<position_type> position;
};

// template<typename SequenceType>
// struct iterator_p
// {
//     using sequence_type = SequenceType;
//     using position_type = position_t<sequence_type>;
//     using iterator_category = std::input_iterator_tag;
//     using value_type        = value_type_t<SequenceType>;
//     using difference_type   = std::ptrdiff_t;
//     using pointer           = value_type*;
//     using reference         = value_type&;

//     iterator_p& operator++() { p = std::make_unique<std::optional<position_type>>(next(sequence, p->value())); return *this; }
//     bool operator!=(const iterator_p& other) const { return p->has_value() != other.p->has_value(); }
//     auto operator*() const { return obtain(sequence, p->value()); }
//     sequence_type sequence;
//     std::unique_ptr<std::optional<position_type>> p;
// };

template<typename SequenceType>
auto begin(SequenceType sequence)
{
    //using ops = std::optional<position_t<SequenceType>>;

    //if constexpr (std::is_copy_assignable<position_t<SequenceType>>::value)
        return iterator<SequenceType>{sequence, start(sequence)};
    //else
    //    return iterator_p<SequenceType>{sequence, std::make_unique<ops>(start(sequence))};
}

template<typename SequenceType>
auto end(SequenceType sequence)
{
    //using ops = std::optional<position_t<SequenceType>>;

    //if constexpr (std::is_copy_assignable<position_t<SequenceType>>::value)
        return iterator<SequenceType>{sequence, {}};
    //else
    //    return iterator_p<SequenceType>{sequence, std::make_unique<ops>(ops{})};
}

template<typename SequenceType>
bool empty(const SequenceType& sequence)
{
    return ! start(sequence).has_value();
}




/**
 * @brief      An infinite sequence containing the values {f(x), f(f(x)), ...}
 *             for a function f: ValueType -> ValueType and an initial ValueType
 *             x.
 *
 * @tparam     ValueType     The type of the sequence values
 * @tparam     FunctionType  The type of the generating function
 */
template<typename ValueType, typename FunctionType>
struct generator_sequence_t
{
    ValueType start;
    FunctionType function;
};

template<typename ValueType, typename FunctionType, typename SequenceType = generator_sequence_t<ValueType, FunctionType>>
std::optional<ValueType> start(generator_sequence_t<ValueType, FunctionType> sequence)
{
    return sequence.start;
}

template<typename ValueType, typename FunctionType, typename SequenceType = generator_sequence_t<ValueType, FunctionType>>
std::optional<ValueType> next(generator_sequence_t<ValueType, FunctionType> sequence, position_t<SequenceType> position)
{
    return sequence.function(position);
}

template<typename ValueType, typename FunctionType, typename SequenceType = generator_sequence_t<ValueType, FunctionType>>
ValueType obtain(generator_sequence_t<ValueType, FunctionType> sequence, position_t<SequenceType> position)
{
    return position;
}

template<typename StartType, typename FunctionType>
auto generate(StartType start, FunctionType function)
{
    using value_type = std::invoke_result_t<FunctionType, StartType>;
    return generator_sequence_t<value_type, FunctionType>{start, function};
}

inline auto generate(long start=0)
{
    return generate(start, [] (long i) { return i + 1; });
}




/**
 * @brief      A sequence containing the values {f(x)} for each x in a source
 *             sequence X, where f is a mapping to a value type, that may be
 *             different from that of X.
 *
 * @tparam     SequenceType  The type of the source sequence X
 * @tparam     FunctionType  The type of the mapping function f
 */
template<typename SequenceType, typename FunctionType>
struct mapped_sequence_t
{
    SequenceType sequence;
    FunctionType mapping;
};

template<typename SequenceType, typename FunctionType>
auto start(mapped_sequence_t<SequenceType, FunctionType> sequence)
{
    return start(sequence.sequence);
}

template<typename SequenceType, typename FunctionType>
auto next(mapped_sequence_t<SequenceType, FunctionType> sequence, position_t<SequenceType> position)
{
    return next(sequence.sequence, position);
}

template<typename SequenceType, typename FunctionType>
auto obtain(mapped_sequence_t<SequenceType, FunctionType> sequence, position_t<SequenceType> position)
{
    return sequence.mapping(obtain(sequence.sequence, position));
}

template<typename SequenceType, typename FunctionType>
auto map(SequenceType sequence, FunctionType function)
{
    return mapped_sequence_t<SequenceType, FunctionType>{sequence, function};
}




/**
 * @brief      A sequence containing those values of a source sequence X for
 *             which a predicate f evaluates to false.
 *
 * @tparam     SequenceType  The type of the source sequence X
 * @tparam     FunctionType  The type of the predicate function f
 */
template<typename SequenceType, typename FunctionType>
struct filtered_sequence_t
{
    SequenceType sequence;
    FunctionType predicate;
};

template<typename SequenceType, typename FunctionType>
auto __not_removed(filtered_sequence_t<SequenceType, FunctionType> sequence, std::optional<position_t<SequenceType>> p)
-> std::optional<position_t<SequenceType>>
{
    if (! p.has_value())
        return {};
    if (! p.has_value() || ! sequence.predicate(obtain(sequence.sequence, p.value())))
        return p;
    return __not_removed(sequence, next(sequence.sequence, p.value()));
}

template<typename SequenceType, typename FunctionType>
auto start(filtered_sequence_t<SequenceType, FunctionType> sequence)
{
    return __not_removed(sequence, start(sequence.sequence));
}

template<typename SequenceType, typename FunctionType>
auto next(filtered_sequence_t<SequenceType, FunctionType> sequence, position_t<SequenceType> position)
{
    return __not_removed(sequence, next(sequence.sequence, position));
}

template<typename SequenceType, typename FunctionType>
auto obtain(filtered_sequence_t<SequenceType, FunctionType> sequence, position_t<SequenceType> position)
{
    return obtain(sequence.sequence, position);
}

template<typename SequenceType, typename FunctionType>
auto remove_if(SequenceType sequence, FunctionType predicate)
{
    return filtered_sequence_t<SequenceType, FunctionType>{sequence, predicate};
}




/**
 * @brief      A sequence containing those values of a source sequence X
 *             preceding the first element for which a predicate function f
 *             evaluates to false.
 *
 * @tparam     SequenceType  The type of the source sequence X
 * @tparam     FunctionType  The type of the predicate function f
 */
template<typename SequenceType, typename FunctionType>
struct truncated_sequence_t
{
    SequenceType sequence;
    FunctionType predicate;
};

template<typename SequenceType, typename FunctionType>
auto start(truncated_sequence_t<SequenceType, FunctionType> sequence)
-> std::optional<position_t<SequenceType>>
{
    if (auto p = start(sequence.sequence); p.has_value() && sequence.predicate(obtain(sequence, p.value())))
        return p;
    return {};
}

template<typename SequenceType, typename FunctionType>
auto next(truncated_sequence_t<SequenceType, FunctionType> sequence, position_t<truncated_sequence_t<SequenceType, FunctionType>> position)
-> std::optional<position_t<SequenceType>>
{
    if (auto p = next(sequence.sequence, position); p.has_value() && sequence.predicate(obtain(sequence, p.value())))
        return p;
    return {};
}

template<typename SequenceType, typename FunctionType>
auto obtain(truncated_sequence_t<SequenceType, FunctionType> sequence, position_t<truncated_sequence_t<SequenceType, FunctionType>> position)
{
    return obtain(sequence.sequence, position);
}

template<typename SequenceType, typename FunctionType>
auto take_while(SequenceType sequence, FunctionType predicate)
{
    return truncated_sequence_t<SequenceType, FunctionType>{sequence, predicate};
}




/**
 * @brief      A sequence that begins at a specified point in a source sequence.
 *
 * @tparam     SequenceType  The type of the source sequence
 */
template<typename SequenceType>
struct advanced_sequence_t
{
    SequenceType sequence;
    std::optional<position_t<SequenceType>> start;
};

template<typename SequenceType>
auto start(advanced_sequence_t<SequenceType> sequence)
{
    return sequence.start;
}

template<typename SequenceType>
auto next(advanced_sequence_t<SequenceType> sequence, position_t<advanced_sequence_t<SequenceType>> position)
{
    return next(sequence.sequence, position);
}

template<typename SequenceType>
auto obtain(advanced_sequence_t<SequenceType> sequence, position_t<advanced_sequence_t<SequenceType>> position)
{
    return obtain(sequence.sequence, position);
}

template<typename SequenceType>
auto advance(SequenceType sequence, std::optional<position_t<SequenceType>> position)
{
    return advanced_sequence_t<SequenceType>{sequence, position};
}




/**
 * @brief      A sequence of the positions in a source sequence.
 *
 * @tparam     SequenceType  The type of the source sequence
 */
template<typename SequenceType>
struct measured_sequence_t
{
    SequenceType sequence;
};

template<typename SequenceType>
auto start(measured_sequence_t<SequenceType> sequence)
{
    return start(sequence.sequence);
}

template<typename SequenceType>
auto next(measured_sequence_t<SequenceType> sequence, position_t<measured_sequence_t<SequenceType>> position)
{
    return next(sequence.sequence, position);
}

template<typename SequenceType>
auto obtain(measured_sequence_t<SequenceType> sequence, position_t<measured_sequence_t<SequenceType>> position)
{
    return position;
}

template<typename SequenceType>
auto measure(SequenceType sequence)
{
    return measured_sequence_t<SequenceType>{sequence};
}




/**
 * @brief      A sequence made by concatenating two source sequences A and B.
 *             The only requirement is that a common type exists for the value
 *             types of A and B.
 *
 * @tparam     SequenceType1  The type of the source sequence A
 * @tparam     SequenceType2  The type of the source sequence B
 */
template<typename SequenceType1, typename SequenceType2>
struct chained_sequence_t
{
    std::pair<SequenceType1, SequenceType2> sequences;
};

template<typename SequenceType1, typename SequenceType2>
auto start(chained_sequence_t<SequenceType1, SequenceType2> sequence)
{
    return detail::optional_pair_of_either(start(sequence.sequences.first), start(sequence.sequences.second));
}

template<typename SequenceType1, typename SequenceType2>
auto next(chained_sequence_t<SequenceType1, SequenceType2> sequence, position_t<chained_sequence_t<SequenceType1, SequenceType2>> position)
{
    return detail::optional_pair_of_either(
          position.first.has_value() ? next(sequence.sequences.first,  position.first .value()) : position.first,
        ! position.first.has_value() ? next(sequence.sequences.second, position.second.value()) : position.second);
}

template<typename SequenceType1, typename SequenceType2>
auto obtain(chained_sequence_t<SequenceType1, SequenceType2> sequence, position_t<chained_sequence_t<SequenceType1, SequenceType2>> position)
{
    using value_type = typename std::common_type<value_type_t<SequenceType1>, value_type_t<SequenceType2>>::type;

    return position.first.has_value()
    ? value_type(obtain(sequence.sequences.first,  position.first .value()))
    : value_type(obtain(sequence.sequences.second, position.second.value()));
}

template<typename SequenceType1, typename SequenceType2>
auto chain(SequenceType1 sequence1, SequenceType2 sequence2)
{
    return chained_sequence_t<SequenceType1, SequenceType2>{std::pair(sequence1, sequence2)};
}

template<typename SequenceType1, typename SequenceType2, typename... SequenceTypes>
auto chain(SequenceType1 sequence1, SequenceType2 sequence2, SequenceTypes... sequences)
{
    return chain(chain(sequence1, sequence2), sequences...);
}




/**
 * @brief      A sequence containing the tuples {{x, y, ...}} for values (x, y,
 *             ...) in the source sequences (X, Y, ...). The zipped sequence
 *             terminates when any of the source sequences terminates.
 *
 * @tparam     SequenceTypes  The types of the source sequences X, Y, ...
 */
template<typename... SequenceTypes>
struct zipped_sequence_t
{
    std::tuple<SequenceTypes...> sequences;
};

template<typename... SequenceTypes>
auto start(zipped_sequence_t<SequenceTypes...> sequence)
{
    using namespace detail;
    return optional_tuple(map(sequence.sequences, [] (auto s) { return start(s); }));
}

template<typename... SequenceTypes>
auto next(zipped_sequence_t<SequenceTypes...> sequence, position_t<zipped_sequence_t<SequenceTypes...>> positions)
{
    using namespace detail;
    return optional_tuple(map(zip_tuples(sequence.sequences, positions), [] (auto sq) {
        return next(std::get<0>(sq), std::get<1>(sq));
    }));
}

template<typename... SequenceTypes>
auto obtain(zipped_sequence_t<SequenceTypes...> sequence, position_t<zipped_sequence_t<SequenceTypes...>> positions)
{
    using namespace detail;
    return map(zip_tuples(sequence.sequences, positions), [] (auto sq) {
        return obtain(std::get<0>(sq), std::get<1>(sq));
    });
}

template<typename SequenceType1, typename SequenceType2>
auto zip(SequenceType1 sequence1, SequenceType2 sequence2)
{
    return map(zipped_sequence_t<SequenceType1, SequenceType2>{{sequence1, sequence2}}, [] (auto t) {
        return std::pair(std::get<0>(t), std::get<1>(t));
    });
}

template<typename... SequenceTypes>
auto zip(SequenceTypes... sequences)
{
    return zipped_sequence_t<SequenceTypes...>{{sequences...}};
}




/**
 * @brief      A sequence made by flattening a sequence of sequences.
 *
 * @tparam     SequenceType  The type of the outer sequence.
 */
template<typename SequenceType>
struct flattened_sequence_t
{
    SequenceType sequence;
};

template<typename SequenceType>
auto start(flattened_sequence_t<SequenceType> sequence)
-> std::optional<std::pair<position_t<SequenceType>, position_t<value_type_t<SequenceType>>>>
{
    if (auto m = start(sequence.sequence); m.has_value())
        return std::pair(m.value(), start(obtain(sequence.sequence, m.value())).value());
    return {};
}

template<typename SequenceType>
auto next(flattened_sequence_t<SequenceType> sequence, position_t<flattened_sequence_t<SequenceType>> position)
-> std::optional<std::pair<position_t<SequenceType>, position_t<value_type_t<SequenceType>>>>
{
    auto outer_sequence = sequence.sequence;
    auto inner_sequence = obtain(outer_sequence, position.first);

    if (auto n = next(inner_sequence, position.second); n.has_value())
        return std::pair(position.first, n.value());

    if (auto m = next(outer_sequence, position.first); m.has_value())
        return std::pair(m.value(), start(obtain(outer_sequence, m.value())).value());

    return {};
}

template<typename SequenceType>
auto obtain(flattened_sequence_t<SequenceType> sequence, position_t<flattened_sequence_t<SequenceType>> position)
{
    return obtain(obtain(sequence.sequence, position.first), position.second);
}

template<typename SequenceType>
auto flat(SequenceType sequence)
{
    auto non_empty = remove_if(sequence, [] (const auto& seq) { return empty(seq); });
    return flattened_sequence_t<decltype(non_empty)>{std::move(non_empty)};
}




/**
 * @brief      A sequence implementing a scan operation. Follows the spec for
 *             Haskell's scanl:
 *
 *             scan(A, b, f) = {b, f(b, a0), f(f(b, a0), a1), ...}
 *
 * @tparam     A     The type of sequence being scanned
 * @tparam     B     The type of the starting value (and return type of the
 *                   reducer)
 * @tparam     C     The type of the reducer function
 */
template<typename A, typename B, typename C>
struct scanned_sequence_t
{
    A sequence;
    B start;
    C reducer;
};

template<typename A, typename B, typename C>
auto start(scanned_sequence_t<A, B, C> sequence)
{
    return std::optional(std::pair(sequence.start, start(sequence.sequence)));
}

template<typename A, typename B, typename C>
auto next(scanned_sequence_t<A, B, C> sequence, position_t<scanned_sequence_t<A, B, C>> position)
-> std::optional<position_t<scanned_sequence_t<A, B, C>>>
{
    if (auto q = position.second; q.has_value())
    {
        auto b = obtain(sequence, position);
        auto a = obtain(sequence.sequence, q.value());
        return std::pair(sequence.reducer(b, a), next(sequence.sequence, q.value()));
    }
    return {};
}

template<typename A, typename B, typename C>
auto obtain(scanned_sequence_t<A, B, C> sequence, position_t<scanned_sequence_t<A, B, C>> position)
{
    return position.first;
}

template<typename A, typename B, typename C>
auto scan(A sequence, B start, C reducer)
{
    static_assert(std::is_same_v<std::invoke_result_t<C, B, value_type_t<A>>, B>,
        "The reducer of a scan operation over an event sequence [a] "
        "with start value type b must be (b, a) -> b");
    return scanned_sequence_t<A, B, C>{sequence, start, reducer};
}




/**
 * @brief      A sequence that adapts any data structure implementing the
 *             sequence protocol to a proper sequence. This technique trivially
 *             wraps this foreign sequence's start, next, and obtain methods,
 *             such that the structure is effectively imported into the seq
 *             namespace. Doing this enables argument-dependent-lookup (ADL) to
 *             find the sequence operators without prefixing the seq namespace,
 *             most importantly the begin and end begin functions, making the
 *             foreign struct iterable via range-based for-loop.
 *
 * @tparam     SequenceType  The type of the data structure implementing the
 *                           sequence protocol
 */
template<typename SequenceType>
struct foreign_sequence_t
{
    SequenceType sequence;
};

template<typename SequenceType>
auto start(foreign_sequence_t<SequenceType> sequence)
{
    return start(sequence.sequence);
}

template<typename SequenceType>
auto next(foreign_sequence_t<SequenceType> sequence, position_t<foreign_sequence_t<SequenceType>> position)
{
    return next(sequence.sequence, position);
}

template<typename SequenceType>
auto obtain(foreign_sequence_t<SequenceType> sequence, position_t<foreign_sequence_t<SequenceType>> position)
{
    return obtain(sequence.sequence, position);
}

template<typename SequenceType>
auto adapt(SequenceType sequence)
{
    return foreign_sequence_t<SequenceType>{sequence};
}




//=============================================================================
template<typename ContainerType>
struct container_sequence_t
{
    std::shared_ptr<ContainerType> container;
};

template<typename ContainerType>
auto start(container_sequence_t<ContainerType> sequence)
{
    const auto& c = *sequence.container;
    return std::empty(c)
    ? std::optional<decltype(std::begin(c))>{}
    : std::optional<decltype(std::begin(c))>{std::begin(c)};
}

template<typename ContainerType>
auto next(container_sequence_t<ContainerType> sequence, position_t<container_sequence_t<ContainerType>> position)
{
    ++position;
    const auto& c = *sequence.container;
    return position == std::end(c)
    ? std::optional<position_t<container_sequence_t<ContainerType>>>{}
    : std::optional<position_t<container_sequence_t<ContainerType>>>{position};
}

template<typename ContainerType>
auto obtain(container_sequence_t<ContainerType> sequence, position_t<container_sequence_t<ContainerType>> position)
{
    return *position;
}

template<typename ContainerType>
auto view(ContainerType&& container)
{
    return container_sequence_t<ContainerType>{std::make_shared<ContainerType>(std::move(container))};
}




//=============================================================================
template<typename ValueType>
struct dynamic_sequence_t
{
    struct base
    {
        virtual std::optional<std::any> next(std::any) const = 0;
        virtual std::optional<std::any> start() const = 0;
        virtual ValueType obtain(std::any) const = 0;
    };
    std::shared_ptr<base> impl;
};

template<typename ValueType>
auto start(dynamic_sequence_t<ValueType> sequence)
{
    return sequence.impl->start();
}

template<typename ValueType>
auto next(dynamic_sequence_t<ValueType> sequence, std::any position)
{
    return sequence.impl->next(position);
}

template<typename ValueType>
ValueType obtain(dynamic_sequence_t<ValueType> sequence, std::any position)
{
    return sequence.impl->obtain(position);
}

template<typename SequenceType>
auto to_dynamic(SequenceType sequence)
{
    using value_type = value_type_t<SequenceType>;
    using position_type = position_t<SequenceType>;

    struct dynamic_sequence_impl : dynamic_sequence_t<value_type>::base
    {
        dynamic_sequence_impl(SequenceType sequence) : sequence(sequence) {}
        std::optional<std::any> next(std::any position) const override { return detail::optional_any(seq::next(sequence, std::any_cast<position_type>(position))); }
        std::optional<std::any> start()                 const override { return detail::optional_any(seq::start(sequence)); }
        value_type obtain(std::any position)            const override { return seq::obtain(sequence, std::any_cast<position_type>(position)); }
        SequenceType sequence;
    };
    return dynamic_sequence_t<value_type>{std::make_shared<dynamic_sequence_impl>(sequence)};
}




//=============================================================================
template<typename SequenceType>
auto front(SequenceType sequence)
{
    return obtain(sequence, start(sequence).value());
}

template<typename SequenceType>
auto back(SequenceType sequence)
{
    auto position = start(sequence).value();

    while (true)
        if (auto p = next(sequence, position); ! p.has_value())
            break;
        else
            position = p.value();

    return obtain(sequence, position);
}

template<typename SequenceType>
auto enumerate(SequenceType sequence)
{
    return zip(generate(), sequence);
}

template<typename SequenceType>
auto take(SequenceType sequence, unsigned long count)
{
    return map(take_while(enumerate(sequence),
        [count] (auto iv) { return (unsigned long)iv.first < count; }),
        []      (auto iv) { return iv.second; });
}

template<typename SequenceType>
auto drop(SequenceType sequence, unsigned long count)
{
    return advance(sequence, std::optional(back(take(measure(sequence), count + 1))));
}

template<typename SequenceType>
auto chunk(SequenceType sequence, unsigned long chunk_size)
{
    auto second = [] (auto s) { return map(s, [] (auto iv) { return iv.second; }); };
    auto not_dividing = [c=chunk_size] (auto iv) { return iv.first % c != 0; };
    auto advance_from = [c=chunk_size, s=sequence] (auto p) { return take(advance(s, std::optional(p)), c); };
    return map(second(remove_if(enumerate(measure(sequence)), not_dividing)), advance_from);
}

inline auto range(unsigned long count)
{
    return take(generate(0), count);
}

inline auto range(long start, long final)
{
    if (final < start)
        throw std::invalid_argument("seq::range (final must be >= start)");
    return take(generate(start), final - start);
}

template<typename... ValueType>
auto from(ValueType... values)
{
    return view(std::array{values...});
}

template<template<typename...> typename ContainerType, typename SequenceType>
auto to(SequenceType sequence)
{
    return ContainerType<value_type_t<SequenceType>>(begin(sequence), end(sequence));
}

template<typename ValueType>
auto always(ValueType value)
{
    return generate(value, [] (auto v) { return v; });
}

template<typename ValueType>
auto just(ValueType value)
{
    return take(always(value), 1);
}

template<typename SequenceType>
auto cycle(SequenceType sequence)
{
    return flat(always(sequence));
}

template<typename SequenceType>
auto repeat(SequenceType sequence, unsigned long count)
{
    return flat(take(always(sequence), count));
}

template<typename SequenceType, typename FunctionType>
auto flat_map(SequenceType sequence, FunctionType function)
{
    return flat(map(sequence, function));
}

template<typename ValueType>
auto yield_if(ValueType value, bool condition)
{
    return take(just(value), condition);
}




/**
 * @brief      Return windowed pairs of the sequence A:
 *
 *             {(a0, a1), (a1, a2), ...}
 *
 * @param[in]  sequence      The source sequence A
 *
 * @tparam     SequenceType  The type of the source sequence A
 *
 * @return     A sequence of pairs
 */
template<typename SequenceType>
auto window(SequenceType sequence)
{
    auto p0 = start(sequence).value();
    auto p1 = next(sequence, p0).value();
    auto a0 = obtain(sequence, p0);
    auto a1 = obtain(sequence, p1);
    auto A2 = advance(sequence, next(sequence, p1));

    return scan(A2, std::pair(a0, a1), [] (auto last_two, auto next_value)
    {
        return std::pair(last_two.second, next_value);
    });
}




//=============================================================================
template<typename SequenceType, typename OperatorType>
auto operator|(SequenceType&& sequence, OperatorType&& op)
{
    return std::forward<OperatorType>(op)(std::forward<SequenceType>(sequence));
}

template<typename T, typename F> auto scan(T b, F f) { return [=] (auto s) { return scan(s, b, f); }; }
template<typename F> auto map       (F f)  { return [f] (auto s) { return map       (s, f); }; }
template<typename F> auto take_while(F f)  { return [f] (auto s) { return take_while(s, f); }; }
template<typename F> auto remove_if (F f)  { return [f] (auto s) { return remove_if (s, f); }; }
template<typename T> auto pair_with (T t)  { return [t] (auto s) { return zip       (s, t); }; }

inline auto flat  ()                { return [ ] (auto s) { return flat  (s); }; }
inline auto cycle ()                { return [ ] (auto s) { return cycle (s); }; }
inline auto window()                { return [ ] (auto s) { return window(s); }; }
inline auto take  (unsigned long c) { return [c] (auto s) { return take  (s, c); }; }
inline auto drop  (unsigned long c) { return [c] (auto s) { return drop  (s, c); }; }
inline auto chunk (unsigned long c) { return [c] (auto s) { return chunk (s, c); }; }
inline auto repeat(unsigned long c) { return [c] (auto s) { return repeat(s, c); }; }

} // namespace seq




//=============================================================================
#ifdef DO_UNIT_TESTS
#include <vector>
#include <string>
#include "core_unit_test.hpp"




//=============================================================================
inline void test_sequence()
{
    require((seq::to<std::basic_string>(seq::view(std::string("tree"))) == std::string("tree")));
    require((seq::to<std::vector>(seq::from(1, 2, 3)) == std::vector{1, 2, 3}));
    require((seq::to<std::vector>(take(cycle(seq::from(1, 2)), 4)) == std::vector{1, 2, 1, 2}));
    require((seq::to<std::vector>(flat(map(seq::range(4), [] (auto i) { return seq::range(i); }))) == std::vector<long>{0, 0, 1, 0, 1, 2}));
    require((seq::to<std::vector>(drop(seq::range(5, 10), 2)) == std::vector{7l, 8l, 9l}));
    require((seq::to<std::vector>(back(chunk(seq::range(6), 3))) == std::vector{3l, 4l, 5l}));
    require((seq::to<std::vector>(scan(seq::from(4, 2, 4), 64.0, std::divides<>())) == std::vector{64., 16., 8., 2.}));
    require((seq::to<std::vector>(window(seq::from(1, 2, 3, 4))) == std::vector{std::pair(1, 2), std::pair(2, 3), std::pair(3, 4)}));
}

#endif // DO_UNIT_TESTS
