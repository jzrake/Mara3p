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
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include "core_numeric_array.hpp"




namespace mesh
{




/**
 * @brief      A struct that identifies a node's global position in the tree:
 *             its level, and its coordinates with respect to the origin at its
 *             level.
 *
 * @tparam     Rank  The rank of the tree to be indexed (having ratio 2^Rank)
 */
template<unsigned long Rank>
struct block_index_t
{
    unsigned long level = 0;
    numeric::array_t<unsigned long, Rank> coordinates = {};
};




/**
 * @brief      Determine if this is a valid index for a binary tree (whether it
 *             is in-bounds on its level).
 *
 * @param[in]  i     The index
 *
 * @tparam     R     The rank
 *
 * @return     True or false
 */
template<unsigned long R>
bool valid_binary_index(const block_index_t<R>& i)
{
    return all(map(i.coordinates, [level=i.level] (auto j) { return j < unsigned(1 << level); }));
}




/**
 * @brief      Return true if this index is at level 0 (in which case an
 *             exception is thrown if the coordinates were not all zero).
 *
 * @param[in]  i     The index
 *
 * @tparam     R     The rank
 *
 * @return     True or false
 */
template<unsigned long R>
bool check_root(const block_index_t<R>& i)
{
    if (i.level == 0)
    {
        if (any(i.coordinates))
        {
            throw std::invalid_argument("mesh::block_index_t (invalid tree index)");
        }
        return true;
    }
    return false;
}




/**
 * @brief      Transform this index, so that it points to the same node as it
 *             does now, but as seen by the child of this node which is in the
 *             direction of the target node.
 *
 * @param[in]  i     The index
 *
 * @tparam     R     The rank
 *
 * @return     The index with level - 1 and the coordinates offset according to
 *             the orthant value
 */
template<unsigned long R>
block_index_t<R> advance_level(const block_index_t<R>& i)
{
    return {i.level - 1, i.coordinates - orthant(i) * (1 << (i.level - 1))};
}




/**
 * @brief      Return this index as it is relative to the parent block.
 *
 * @param[in]  i     The index
 *
 * @tparam     R     The rank
 *
 * @return     The index as seen by the parent.
 */
template<unsigned long R>
block_index_t<R> relative_to_parent(const block_index_t<R>& i)
{
    return {1, map(i.coordinates, [] (auto c) { return c % 2; })};
}




/**
 * @brief      Return the tree index at level - 1 which contains this one as a
 *             child.
 *
 * @param[in]  i     The index
 *
 * @tparam     R     The rank
 *
 * @return     The parent index
 */
template<unsigned long R>
block_index_t<R> parent_index(const block_index_t<R>& i)
{
    return {i.level - 1, i.coordinates / 2};
}




/**
 * @brief      Return a sequence of tree indexes corresponding to the 2^R
 *             indexes of this node's children.
 *
 * @param[in]  i     The index
 *
 * @tparam     R     The rank
 *
 * @return     A sequence of tree indexes
 */
template<unsigned long R>
numeric::array_t<block_index_t<R>, 1 << R> child_indexes(const block_index_t<R>& i)
{
    return map(numeric::range<1 << R>(), [&i] (auto j) -> block_index_t<R>
    {
        return {i.level + 1, i.coordinates * 2 + numeric::binary_repr<R>(j)};
    });
}




/**
 * @brief      Helper function to give the number of children per node at the
 *             given Rank
 *
 * @param[in]  <unnamed>  { parameter_description }
 *
 * @tparam     R          { description }
 *
 * @return     An integer
 */
template<unsigned long R>
std::size_t child_count(const block_index_t<R>&)
{
    return 1 << R;
}




/**
 * @brief      Return the linear index (between 0 and 2^Rank, non-inclusive) of
 *             this index in its parent.
 *
 * @param[in]  i     { parameter_description }
 *
 * @tparam     R     { description }
 *
 * @return     The sibling index
 */
template<unsigned long R>
std::size_t sibling_index(const block_index_t<R>& i)
{
    return to_integral(orthant(relative_to_parent(i)));
}




/**
 * @brief      Determines whether the specified i is last child.
 *
 * @param[in]  i     { parameter_description }
 *
 * @tparam     R     { description }
 *
 * @return     True if the specified i is last child, False otherwise.
 */
template<unsigned long R>
bool is_last_child(const block_index_t<R>& i)
{
    return i.level == 0 || sibling_index(i) == child_count(i) - 1;
}




/**
 * @brief      Return the next sibling's index. Throws std::out_of_range if
 *             next sibling does not exist
 *
 * @return     An index
 */
template<unsigned long R>
block_index_t<R> next_sibling(const block_index_t<R>& i)
{
    if (is_last_child(i))
    {
        throw std::out_of_range("mesh::next_sibling (is the last child)");
    }
    return block_index_t<R>{i.level, parent_index(i).coordinates * 2 + numeric::binary_repr<R>(sibling_index(i) + 1)};
}




/**
 * @brief      Return a new index with the same coordinates but a different
 *             level.
 *
 * @param[in]  i          The index
 * @param[in]  new_level  The level of the returned index
 *
 * @tparam     R          The rank
 *
 * @return     A new index
 */
template<unsigned long R>
block_index_t<R> with_level(const block_index_t<R>& i, unsigned long new_level)
{
    return {new_level, i.coordinates};
}




/**
 * @brief      Return a new index with the same coordinates but a different
 *             level.
 *
 * @param[in]  i                The index
 * @param[in]  new_coordinates  The coordinates of the returned index
 *
 * @tparam     R                The rank
 *
 * @return     A new index
 */
template<unsigned long R>
block_index_t<R> with_coordinates(const block_index_t<R>& i, numeric::array_t<unsigned long, R> new_coordinates)
{
    return {i.level, new_coordinates};
}




/**
 * @brief      Return the orthant (ray, quadrant, octant) of this index.
 *
 * @param[in]  i     The index
 *
 * @tparam     R     The rank
 *
 * @return     A sequence of bool's
 *
 * @note       https://en.wikipedia.org/wiki/Orthant
 */
template<unsigned long R>
numeric::array_t<bool, R> orthant(const block_index_t<R> &i)
{
    return map(i.coordinates, [level=i.level] (auto x) -> bool { return x / (1 << (level - 1)); });
}




/**
 * @brief      Return an index shifted up on the given axis. The result index is
 *             wrapped between 0 and the maximum value at this level.
 *
 * @param[in]  i     The index
 * @param[in]  a     The axis to shift on
 *
 * @tparam     R     The rank
 *
 * @return     A new index, on the same level
 */
template<unsigned long R> block_index_t<R> next_on(const block_index_t<R>& i, std::size_t a) { return {i.level, update(i.coordinates, a, [L=i.level] (auto j) { return (j + (1 << L) + 1) % (1 << L); })}; }
template<unsigned long R> block_index_t<R> prev_on(const block_index_t<R>& i, std::size_t a) { return {i.level, update(i.coordinates, a, [L=i.level] (auto j) { return (j + (1 << L) - 1) % (1 << L); })}; }




//=============================================================================
template<unsigned long R> bool operator==(const block_index_t<R>& a, const block_index_t<R>& b) { return a.level == b.level && a.coordinates == b.coordinates; }
template<unsigned long R> bool operator!=(const block_index_t<R>& a, const block_index_t<R>& b) { return a.level != b.level || a.coordinates != b.coordinates; }
template<unsigned long R> bool operator< (const block_index_t<R>& a, const block_index_t<R>& b) { return a.level != b.level ? a.level < b.level : a.coordinates.impl < b.coordinates.impl; }
template<unsigned long R> bool operator> (const block_index_t<R>& a, const block_index_t<R>& b) { return a.level != b.level ? a.level > b.level : a.coordinates.impl > b.coordinates.impl; }





/**
 * @brief      Return a stringified tree index that can be re-parsed by the
 *             read_block_index method.
 *
 * @param[in]  index  The tree index
 *
 * @tparam     Rank   The rank of the tree
 *
 * @return     A string, formatted like "level:i-j-k"
 */
template<unsigned long Rank>
std::string format_block_index(block_index_t<Rank> index)
{
    std::stringstream ss;
    ss << index.level;

    for (std::size_t i = 0; i < Rank; ++i)
    {
        ss
        << (i == 0 ? ':' : '-')
        << std::setfill('0')
        << std::setw(1 + std::log10(1 << index.level))
        << index.coordinates[i];
    }
    return ss.str();
}




/**
 * @brief      Return a tree index from a string formatted with
 *             format_block_index.
 *
 * @param[in]  str   The string representation of the index
 *
 * @tparam     Rank  The rank of the tree
 *
 * @return     The index
 */
template<unsigned long Rank>
block_index_t<Rank> read_block_index(std::string str)
{
    if (std::count(str.begin(), str.end(), '-') != Rank - 1)
    {
        throw std::invalid_argument("mesh::read_block_index (string has wrong rank)");
    }
    auto result = block_index_t<Rank>();

    result.level = std::stoi(str.substr(0, str.find(':')));
    str.erase(0, str.find(':') + 1);

    for (std::size_t i = 0; i < Rank; ++i)
    {
        result.coordinates[i] = std::stoi(str.substr(0, str.find('-')));
        str.erase(0, str.find('-') + 1);
    }
    return result;
}




/**
 * @brief      Return a default-constructed block_index_t
 *
 * @tparam     Rank  The rank of the tree
 *
 * @return     A new tree index
 */
template<unsigned long Rank>
block_index_t<Rank> block_index()
{
    return block_index_t<Rank>{};
}

} // namespace mesh
