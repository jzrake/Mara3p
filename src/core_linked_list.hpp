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
#include <optional>
#include <memory>




//=============================================================================
namespace list
{




//=============================================================================
template<typename ValueType>
struct singly_linked_t
{

    singly_linked_t() {}
    singly_linked_t(const singly_linked_t& other) : first(other.first), rest(other.rest) {}
    singly_linked_t(singly_linked_t&& other) : first(std::move(other.first)), rest(std::move(other.rest)) {}
    singly_linked_t(ValueType first, std::shared_ptr<singly_linked_t> rest) : first(first), rest(std::move(rest)) {}


    /**
     * @brief      Construct a linked list from a pair of iterators
     *
     * @param[in]  start     The begin iterator
     * @param[in]  last      The end iterator
     *
     * @tparam     Iterator  The iterator type
     */
    template<typename Iterator>
    singly_linked_t(Iterator start, Iterator last)
    {
        auto result = singly_linked_t();

        while (start != last)
        {
            result = prepend(result, *start);
            ++start;
        }
        *this = reverse(std::move(result));
    }


    /**
     * @brief      Destroy the linked list
     *
     * @note       This destructor makes sure not to overflow the stack by
     *             recursive calls to ~shared_ptr. See Sec. 8.1.4 of Functional
     *             Programming in C++ by Ivan Čukić.
     */
    ~singly_linked_t()
    {
        while (rest.unique())
        {
            rest = std::move(rest->rest);
        }
    }


    singly_linked_t& operator=(const singly_linked_t& other)
    {
        first = other.first;
        rest = other.rest;
        return *this;
    }


    singly_linked_t& operator=(singly_linked_t&& other)
    {
        first = std::move(other.first);
        rest = std::move(other.rest);
        return *this;
    }


    ValueType first;
    std::shared_ptr<singly_linked_t> rest;
};




/**
 * @brief      Return a single linked list of a single value
 *
 * @param[in]  first      The value
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     A new linked list
 */
template<typename ValueType>
auto just(ValueType first)
{
    return singly_linked_t<ValueType>{first, std::make_shared<singly_linked_t<ValueType>>()};
}




/**
 * @brief      Return an empty list with the given value type
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     An empty linked list
 */
template<typename ValueType>
auto empty()
{
    return singly_linked_t<ValueType>{};
}




/**
 * @brief      Return true if this list is empty
 *
 * @param[in]  list       The list to check for emptiness
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     True or false
 */
template<typename ValueType>
bool empty(const singly_linked_t<ValueType>& list)
{
    return list.rest == nullptr;
}




/**
 * @brief      Return true if this list is empty
 *
 * @param[in]  list       The list to check for emptiness
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     True or false
 */
template<typename ValueType>
std::size_t size(const singly_linked_t<ValueType>& list)
{
    auto count = std::size_t(0);
    auto rest = list.rest;

    while (rest)
    {
        ++count;
        rest = rest->rest;
    }
    return count;
}




/**
 * @brief      Reverse this list; O(N).
 *
 * @param[in]  list       The list to reverse
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     A reversed version of the given list
 */
template<typename ValueType>
auto reverse(singly_linked_t<ValueType> A)
{
    auto B = empty<ValueType>();

    while (! empty(A))
    {
        B = prepend(B, head(A));
        A = tail(A);
    }
    return B;
}




/**
 * @brief      { function_description }
 *
 * @param[in]  args        The arguments
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     { description_of_the_return_value }
 */
template<typename... ValueTypes>
auto from(ValueTypes... args)
{
    using value_type = std::common_type_t<ValueTypes...>;
    auto list = empty<value_type>();

    for (auto a : {value_type(args)...})
    {
        list = prepend(list, a);
    }
    return reverse(list);
}




/**
 * @brief      Return the value at the front of the list; O(1). Throws
 *             out_of_range if this list is empty.
 *
 * @param[in]  list       The list to get the head of
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     The value at the front
 */
template<typename ValueType>
const ValueType& head(const singly_linked_t<ValueType>& list)
{
    if (empty(list))
    {
        throw std::out_of_range("list::singly_linked_t (cannot get the head of an empty list)");
    }
    return list.first;
}




/**
 * @brief      Return the rest of the list, if this list is non-empty; O(1).
 *
 * @param[in]  list       The list
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     The rest of the list
 */
template<typename ValueType>
const singly_linked_t<ValueType>& tail(const singly_linked_t<ValueType>& list)
{
    if (empty(list))
    {
        throw std::out_of_range("list::singly_linked_t (cannot get the tail of an empty list)");
    }
    return *list.rest;
}




/**
 * @brief      Return another list with the given value prepended; O(1).
 *
 * @param[in]  rest       The list to append to
 * @param[in]  first      The value to prepend
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     A new singly linked list
 */
template<typename ValueType>
auto prepend(const singly_linked_t<ValueType>& rest, ValueType first)
{
    return singly_linked_t<ValueType>{first, std::make_shared<singly_linked_t<ValueType>>(rest)};
}




/**
 * @brief      Concatenate two lists; O(N).
 *
 * @param[in]  a          The first list to concatenate
 * @param[in]  b          The second list
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     The elements of this list, followed by the elements the another
 */
template<typename ValueType>
auto concat(const singly_linked_t<ValueType>& a, const singly_linked_t<ValueType>& b)
-> singly_linked_t<ValueType>
{
    return empty(a) ? b : prepend(concat(tail(a), b), head(a));
}




/**
 * @brief      Remove elements from the beginning of a list.
 *
 * @param[in]  list       The list to drop elements from
 * @param[in]  count      The number of elements to drop
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     A list, shortened from the front
 */
template<typename ValueType>
auto drop(singly_linked_t<ValueType> list, std::size_t count)
{
    while (count--)
    {
        if (empty(list))
        {
            throw std::out_of_range("list::drop (cannot drop more elements than list size)");
        }
        list = tail(list);
    }
    return list;
}




/**
 * @brief      { function_description }
 *
 * @param[in]  list       The list
 * @param[in]  count      The count
 *
 * @tparam     ValueType  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType>
auto take(const singly_linked_t<ValueType>& list, std::size_t count)
{
    return reverse(drop(reverse(list), size(list) - count));
}




/**
 * @brief      Sort and merge two already sorted lists.
 *
 * @param[in]  a               The first list
 * @param[in]  b               The second list
 * @param      cmp             A comparator function object
 *
 * @tparam     ValueType       The value type of the list
 * @tparam     ComparatorType  The type of the comparator function
 *
 * @return     A sorted list
 */
template<typename ValueType, typename ComparatorType>
auto sorted_merge(const singly_linked_t<ValueType>& a, const singly_linked_t<ValueType>& b, ComparatorType&& cmp)
{
    if (empty(a)) return b;
    if (empty(b)) return a;

    return cmp(head(a), head(b))
    ? prepend(sorted_merge(tail(a), b, cmp), head(a))
    : prepend(sorted_merge(tail(b), a, cmp), head(b));
}




/**
 * @brief      { function_description }
 *
 * @param[in]  list            The list
 * @param      cmp             The compare
 *
 * @tparam     ValueType       { description }
 * @tparam     ComparatorType  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ComparatorType>
singly_linked_t<ValueType> sort(const singly_linked_t<ValueType>& list, ComparatorType&& cmp)
{
    auto n = size(list);

    if (n == 0 || n == 1)
    {
        return list;
    }
    return sorted_merge(sort(take(list, n / 2), cmp), sort(drop(list, n / 2), cmp), cmp);
}




/**
 * @brief      { function_description }
 *
 * @param[in]  list       The list
 *
 * @tparam     ValueType  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType>
singly_linked_t<ValueType> sort(const singly_linked_t<ValueType>& list)
{
    return sort(list, std::less<>());
}




/**
 * @brief      Create a list of all unique values in a list
 *
 * @return     A sorted list of unique values
 */
template<typename ValueType>
singly_linked_t<ValueType> unique(const singly_linked_t<ValueType>& list)
{
    auto sorted = sort(list);
    auto prev   = head(sorted);
    auto result = just(prev);

    while (! empty(sorted))
    {
        if (prev != head(sorted))
        {
            result = prepend(result, head(sorted));
        }
        prev   = head(sorted);
        sorted = tail(sorted);
    }
    return reverse(result);
}




/**
 * @brief      Determine whether two sequences are equal.
 *
 * @param[in]  a          { parameter_description }
 * @param[in]  b          { parameter_description }
 *
 * @tparam     ValueType  { description }
 *
 * @return     True or false
 */
template<typename ValueType>
bool operator==(const singly_linked_t<ValueType>& a, const singly_linked_t<ValueType>& b)
{
    if (! empty(a) && ! empty(b))
    {
        return head(a) == head(b) && tail(a) == tail(b);
    }
    if (empty(a) && empty(b))
    {
        return true;
    }
    return false;
}




/**
 * @brief      Determine whether two sequences are unequal.
 *
 * @param[in]  a          { parameter_description }
 * @param[in]  b          { parameter_description }
 *
 * @tparam     ValueType  { description }
 *
 * @return     True or false
 */
template<typename ValueType>
bool operator!=(const singly_linked_t<ValueType>& a, const singly_linked_t<ValueType>& b)
{
    return ! operator==(a, b);
}




/**
 * @brief      Functions implementing the sequence protocol. These make it
 *             possible to adapt a linked list to a sequence via seq::adapt.
 *
 * @param[in]  list       The linked list
 *
 * @tparam     ValueType  The value type of the list
 *
 * @return     An optional to the list start, as expected by sequence operators
 */
template<typename ValueType>
auto start(const singly_linked_t<ValueType>& list) -> std::optional<singly_linked_t<ValueType>>
{
    if (! empty(list))
    {
        return list;
    }
    return {};
}

template<typename ValueType>
auto next(const singly_linked_t<ValueType>&, const singly_linked_t<ValueType>& list) -> std::optional<singly_linked_t<ValueType>>
{
    if (! empty(tail(list)))
    {
        return tail(list);
    }
    return {};
}

template<typename ValueType>
auto obtain(const singly_linked_t<ValueType>& a, const singly_linked_t<ValueType>& b)
{
    return b.first;
}

} // namespace list




//=============================================================================
#ifdef DO_UNIT_TESTS
#include <vector>
#include "core_sequence.hpp"
#include "core_unit_test.hpp"




//=============================================================================
inline void test_linked_list()
{
    require(empty(list::empty<int>()));
    require(! empty(list::just(12)));
    require(head(list::just(12)) == 12);
    require(head(list::from(12, 13, 14)) == 12);
    require(head(prepend(list::just(12), 13)) == 13);
    require(size(list::from(12, 13, 14)) == 3);

    auto values = std::vector{1., 2., 3., 4., 5., 6.};
    auto big = 1000;

    require((seq::to<std::vector>(seq::adapt(list::singly_linked_t<double>(values.begin(), values.end()))) == values));
    require(head(seq::to<list::singly_linked_t>(seq::range(10))) == 0);
    require(size(seq::to<list::singly_linked_t>(seq::range(10))) == 10);
    require(size(seq::to<list::singly_linked_t>(seq::range(big))) == std::size_t(big));
    require(size(take(seq::to<list::singly_linked_t>(seq::range(10)), 2)) == std::size_t(2));
    require(head(take(seq::to<list::singly_linked_t>(seq::range(10)), 2)) == std::size_t(0));


    require(head(drop(seq::to<list::singly_linked_t>(seq::range(10)), 2)) == 2);
    require(empty(drop(seq::to<list::singly_linked_t>(seq::range(10)), 10)));
    require_throws(drop(seq::to<list::singly_linked_t>(seq::range(10)), 11));


    require(head(drop(seq::to<list::singly_linked_t>(seq::range(big)), big - 1)) + 1 == big);
    require(head(reverse(seq::to<list::singly_linked_t>(seq::range(big)))) == big - 1);
    require(sorted_merge(list::from(1, 3, 5), list::from(2, 4, 6), std::less<>()) == list::from(1, 2, 3, 4, 5, 6));


    auto A = list::from(1., 2., 3., 4., 5., 6.);
    auto B = concat(take(A, 3), drop(A, 3));
    require(A == B);
    require(seq::to<std::vector>(seq::adapt(sort(list::from(1., 2., 3., 4., 5., 6.), std::less<>()))) == values);

    require(unique(list::from(1, 2, 3)) == list::from(1, 2, 3));
    require(unique(list::from(3, 2, 2, 1, 2, 1, 3, 3)) == list::from(1, 2, 3));
}

#endif // DO_UNIT_TESTS
