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
#include "core_bsp_tree.hpp"
#include "core_numeric_array.hpp"
#include "mesh_block_index.hpp"




//=============================================================================
namespace bsp { // bqo := binary tree, quad-tree, oct-tree




//=============================================================================
template<unsigned long Ratio> struct log2 {};
template<> struct log2<2> { static const unsigned long value = 1; };
template<> struct log2<4> { static const unsigned long value = 2; };
template<> struct log2<8> { static const unsigned long value = 3; };




/**
 * @brief      Return the child of the given (non-leaf) tree node. This function
 *             wraps the bsp method based on linear indexes.
 *
 * @param[in]  tree          The tree (throws out_of_range if not a leaf)
 * @param[in]  orthant       The orthant (ray, quadrant, octant)
 *
 * @tparam     ValueType     The value type of the tree
 * @tparam     ChildrenType  The tree's children provider type
 * @tparam     Ratio         The ratio of the tree (2, 4, or 8)
 *
 * @return     A tree node that is a child of the given node
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto child_at(const tree_t<ValueType, ChildrenType, Ratio>& tree, numeric::array_t<bool, log2<Ratio>::value> orthant)
{
    return child_at(tree, to_integral(orthant));
}




/**
 * @brief      Return a node at an arbitrary depth below the given node. If no
 *             node exists there, throws out_of_range.
 *
 * @param[in]  tree          The tree to index
 * @param[in]  index         The index into the tree
 *
 * @tparam     ValueType     The value type of the tree
 * @tparam     ChildrenType  The tree's children provider type
 * @tparam     Ratio         The ratio of the tree (2, 4, or 8)
 *
 * @return     A tree node that is a descendant of the given node
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto node_at(const tree_t<ValueType, ChildrenType, Ratio>& tree, mesh::block_index_t<log2<Ratio>::value> index)
-> tree_t<ValueType, ChildrenType, Ratio>
{
    return check_root(index) ? tree : node_at(child_at(tree, orthant(index)), advance_level(index));
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree          The tree
 * @param[in]  index         The index
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto value_at(const tree_t<ValueType, ChildrenType, Ratio>& tree, mesh::block_index_t<log2<Ratio>::value> index)
{
    return value(node_at(tree, index));
}




/**
 * @brief      Determine whether a leaf node exists at an arbitrary depth below
 *             the given node.
 *
 * @param[in]  tree          The tree to check
 * @param[in]  index         The index into the tree
 *
 * @tparam     ValueType     The value type of the tree
 * @tparam     ChildrenType  The tree's children provider type
 * @tparam     Ratio         The ratio of the tree
 *
 * @return     True or false
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
bool contains(const tree_t<ValueType, ChildrenType, Ratio>& tree, mesh::block_index_t<log2<Ratio>::value> index)
{
    return has_value(tree) ? check_root(index) : contains(child_at(tree, orthant(index)), advance_level(index));
}




/**
 * @brief      Determine whether there is a node (leaf or otherwise) at the
 *             given index.
 *
 * @param[in]  tree          The tree to check
 * @param[in]  index         The index into the tree
 *
 * @tparam     ValueType     The value type of the tree
 * @tparam     ChildrenType  The tree's children provider type
 * @tparam     Ratio         The ratio of the tree
 *
 * @return     True or false
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
bool contains_node(const tree_t<ValueType, ChildrenType, Ratio>& tree, mesh::block_index_t<log2<Ratio>::value> index)
{
    if (index.level == 0)
    {
        return ! any(index.coordinates);
    }
    if (has_value(tree))
    {
        return false;
    }
    return contains_node(child_at(tree, orthant(index)), advance_level(index));
}




/**
 * @brief      Insert or replace a value at the given index, creating
 *             intermediate nodes as necessary. Throws an exception if a
 *             non-leaf node already exists at the target index.
 *
 * @param[in]  tree       The tree to insert a node to
 * @param[in]  index      The target index
 * @param[in]  value      The value to insert at that index
 *
 * @tparam     ValueType  The tree value type
 * @tparam     Ratio      The tree ratio
 *
 * @return     A new tree
 *
 * @note       This method may generate nodes with default-constructed values at
 *             indexes other than the target index, so the value type must be
 *             default-constructible. Its most likely use is in loading data
 *             into the tree from a file.
 */
template<typename ValueType, unsigned long Ratio>
shared_tree<ValueType, Ratio> insert(shared_tree<ValueType, Ratio> tree, mesh::block_index_t<log2<Ratio>::value> index, ValueType value)
{
    if (index.level == 0)
    {
        if (any(index.coordinates) || ! has_value(tree))
        {
            throw std::out_of_range("bsp::insert (target node already has children)");
        }
        return just<Ratio>(value);
    }
    if (has_value(tree))
    {
        return insert(from(numeric::array_t<ValueType, Ratio>()), index, value);
    }
    return {shared_trees(update(*children(tree).ptr, to_integral(orthant(index)), [next=advance_level(index), value] (auto c)
    {
        return insert(c, next, value);
    }))};
}




/**
 * @brief      Return a tree of indexes representing the topology of the given
 *             tree.
 *
 * @param[in]  tree             The tree node
 * @param[in]  index_in_parent  The starting index (optional)
 *
 * @tparam     ValueType        The tree value type
 * @tparam     ChildrenType     The tree provider type
 * @tparam     Ratio            The tree ratio
 * @tparam     Rank             The rank of the tree: log2(ratio)
 *
 * @return     A tree of indexes
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, unsigned long Rank=log2<Ratio>::value>
auto indexes(const tree_t<ValueType, ChildrenType, Ratio>& tree, mesh::block_index_t<Rank> index_in_parent={})
{
    using provider_type = shared_children_t<mesh::block_index_t<Rank>, Ratio>;

    if (has_value(tree))
    {
        return just<Ratio>(index_in_parent);
    }
    auto A = child_indexes(index_in_parent);
    auto B = children(tree);
    auto C = map(numeric::range<Ratio>(), [A, B] (auto i) { return indexes(B(i), A[i]); });

    return tree_t<mesh::block_index_t<Rank>, provider_type, Ratio>{shared_trees(C)};
}





/**
 * @brief      Like enumerate: zip(indexes(tree), tree).
 *
 * @param[in]  tree          The tree to indexify
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     A tree of two-tuples
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto indexify(const tree_t<ValueType, ChildrenType, Ratio>& tree)
{
    return zip(indexes(tree), tree);
}




/**
 * @brief      Functions implementing the sequence protocol. These make it
 *             possible to adapt a bqo tree to a sequence via seq::adapt.
 *
 * @param[in]  tree          The tree to adapt to a sequence
 * @param[in]  current       The starting index (used internally)
 *
 * @tparam     ValueType     The value type of the tree
 * @tparam     ChildrenType  The tree's children provider type
 *
 * @return     An optional to the tree start, as expected by sequence operators
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, unsigned long Rank=log2<Ratio>::value>
std::optional<mesh::block_index_t<Rank>> start(tree_t<ValueType, ChildrenType, Ratio> tree, mesh::block_index_t<Rank> current={})
{
    return has_value(tree) ? current : start(child_at(tree, 0), child_indexes(current).at(0));
}

template<typename ValueType, typename ChildrenType, unsigned long Ratio, unsigned long Rank=log2<Ratio>::value>
std::optional<mesh::block_index_t<Rank>> next(tree_t<ValueType, ChildrenType, Ratio> tree, mesh::block_index_t<Rank> current)
{
    if (check_root(current))
    {
        return {};
    }

    if (is_last_child(current))
    {
        auto ancestor = parent_index(current);

        while (is_last_child(ancestor))
        {
            if (check_root(ancestor))
            {
                return {};
            }
            ancestor = parent_index(ancestor);
        }
        current = next_sibling(ancestor);

        while (! contains(tree, current))
        {
            current = child_indexes(current).at(0);
        }
        return current;
    }
    current = next_sibling(current);

    while (! contains(tree, current))
    {
        current = child_indexes(current).at(0);
    }
    return current;
}

template<typename ValueType, typename ChildrenType, unsigned long Ratio, unsigned long Rank=log2<Ratio>::value>
auto obtain(tree_t<ValueType, ChildrenType, Ratio> tree, mesh::block_index_t<Rank> index)
{
    return value(node_at(tree, index));
}




//=============================================================================
template<unsigned long Rank>
auto uniform_mesh_tree(unsigned long depth)
{
    auto branch_function = [] (auto i) { return child_indexes(i); };
    auto result = just<1 << Rank>(mesh::block_index_t<Rank>());

    while (depth--)
    {
        result = branch_all(result, branch_function);
    }
    return result;
}




//=============================================================================
template<unsigned long Rank>
inline auto mesh_tree(std::function<bool(mesh::block_index_t<Rank>)> predicate, unsigned long max_depth)
{
    auto branch_function = [] (auto i) { return child_indexes(i); };
    auto result = just<1 << Rank>(mesh::block_index_t<Rank>());

    while (max_depth--)
    {
        result = branch_if(result, branch_function, predicate);
    }
    return result;
}

} // namespace bsp




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"
#include "core_sequence.hpp"




//=============================================================================
inline void test_bqo_tree()
{
    require(orthant(mesh::block_index_t<3>{1, {0, 0, 0}}) == numeric::array(false, false, false));
    require(orthant(mesh::block_index_t<3>{1, {0, 0, 1}}) == numeric::array(false, false, true));
    require(orthant(mesh::block_index_t<3>{1, {0, 1, 0}}) == numeric::array(false, true, false));
    require(orthant(mesh::block_index_t<3>{1, {1, 0, 0}}) == numeric::array(true, false, false));
    require(orthant(mesh::block_index_t<3>{2, {0, 0, 1}}) == numeric::array(false, false, false));
    require(orthant(mesh::block_index_t<3>{2, {0, 1, 0}}) == numeric::array(false, false, false));
    require(orthant(mesh::block_index_t<3>{2, {1, 0, 0}}) == numeric::array(false, false, false));
    require(orthant(mesh::block_index_t<3>{2, {0, 0, 2}}) == numeric::array(false, false, true));
    require(orthant(mesh::block_index_t<3>{2, {0, 2, 0}}) == numeric::array(false, true, false));
    require(orthant(mesh::block_index_t<3>{2, {2, 0, 0}}) == numeric::array(true, false, false));

    require_throws(next_on(mesh::block_index_t<3>(), 3));
    require(  next_on(with_level(mesh::block_index_t<3>(), 3), 0).coordinates == numeric::array(1, 0, 0));
    require(  next_on(with_level(mesh::block_index_t<3>(), 3), 1).coordinates == numeric::array(0, 1, 0));
    require(  next_on(with_level(mesh::block_index_t<3>(), 3), 2).coordinates == numeric::array(0, 0, 1));
    require(  prev_on(with_level(mesh::block_index_t<3>(), 3), 0).coordinates == numeric::array(7, 0, 0));
    require(  prev_on(with_level(mesh::block_index_t<3>(), 3), 1).coordinates == numeric::array(0, 7, 0));
    require(  prev_on(with_level(mesh::block_index_t<3>(), 3), 2).coordinates == numeric::array(0, 0, 7));
    require(  valid_binary_index(with_coordinates(with_level(mesh::block_index_t<3>(), 3), {0, 0, 7})));
    require(! valid_binary_index(with_coordinates(with_level(mesh::block_index_t<3>(), 3), {0, 0, 8})));

    auto branch_function = [] (auto i) { return child_indexes(i); };
    auto tree64 = branch_all(branch_all(bsp::just<8>(mesh::block_index_t<3>()), branch_function), branch_function);
    auto index = with_coordinates(with_level(mesh::block_index_t<3>(), 2), {0, 1, 2});

    require(size(bsp::from(child_indexes(mesh::block_index_t<3>()))) == 8);
    require(size(branch(bsp::just<8>(mesh::block_index_t<3>()), branch_function)) == 8);
    require(size(tree64) == 64);
    require(has_value(bsp::node_at(tree64, index)));
    require(  contains(tree64, index));
    require(! contains(tree64, child_indexes(index)[0]));
    require(! contains(tree64, parent_index(index)));
    require_throws(bsp::node_at(tree64, child_indexes(index)[0]));
}

#endif // DO_UNIT_TESTS
