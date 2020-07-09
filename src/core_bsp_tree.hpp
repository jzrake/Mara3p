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
#include <memory>
#include <variant>
#include <array>
#include "core_numeric_array.hpp"




// bsp := Binary Space Partitioning (facilitates kd-tree, quad-tree, oct-tree, etc.)
namespace bsp {




//=============================================================================
namespace detail {

template<typename FunctionType, typename... Ts>
auto map(std::tuple<Ts...> t, FunctionType fn)
{
    return apply([fn] (auto... ts) { return std::tuple(fn(ts)...); }, t);
}

template<std::size_t Size, typename FunctionType>
auto sum(FunctionType&& f)
{
    auto result = std::invoke_result_t<FunctionType, std::size_t>(0);

    for (std::size_t i = 0; i < Size; ++i)
        result += f(i);

    return result;
}

} // namespace detail




//=============================================================================
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
struct tree_t
{
    std::variant<ValueType, ChildrenType> provider;
};




//=============================================================================
template<typename ValueType, unsigned long Ratio>
struct shared_children_t
{
    static constexpr unsigned long ratio = Ratio;
    using value_type = ValueType;

    auto operator()(std::size_t i) const
    {
        return ptr->operator[](i);
    }
    std::shared_ptr<numeric::array_t<tree_t<value_type, shared_children_t, ratio>, ratio>> ptr;
};

template<typename ChildrenType, typename FunctionType>
struct mapped_children_t
{
    auto operator()(std::size_t i) const
    {
        return map(children(i), function);
    }
    ChildrenType children;
    FunctionType function;
};

template<typename... ChildrenType>
struct zipped_children_t
{
    auto operator()(std::size_t i) const
    {       
        return apply([] (auto... c) { return zip(c...); }, detail::map(children, [i] (auto c) { return c(i); }));
    }
    std::tuple<ChildrenType...> children;
};




//=============================================================================
template<typename ValueType, unsigned long Ratio>
using shared_tree = tree_t<ValueType, shared_children_t<ValueType, Ratio>, Ratio>;

template<unsigned long Ratio, typename ValueType>
auto just(ValueType value)
{
    return shared_tree<ValueType, Ratio>{value};
}

template<typename ValueType, unsigned long Ratio>
auto shared_trees(numeric::array_t<shared_tree<ValueType, Ratio>, Ratio> child_trees)
{
    return shared_children_t<ValueType, Ratio>{
        std::make_shared<decltype(child_trees)>(child_trees)
    };
}

template<typename ValueType, unsigned long Ratio>
auto shared_values(numeric::array_t<ValueType, Ratio> child_values)
{
    return shared_trees(map(child_values, [] (auto v) { return just<Ratio>(v); }));
}

template<typename... ValueTypes>
auto shared_values(ValueTypes... values)
{
    return shared_values(numeric::array(values...));
}

template<typename... ValueTypes>
auto from(ValueTypes... values)
{
    using value_type = std::common_type_t<ValueTypes...>;
    return shared_tree<value_type, sizeof...(ValueTypes)>{shared_values(numeric::array(values...))};
}

template<typename ValueType, std::size_t Ratio>
auto from(numeric::array_t<ValueType, Ratio> values)
{
    return shared_tree<ValueType, Ratio>{shared_values(values)};
}




//=============================================================================
template<typename ValueType, unsigned long Ratio, typename FunctionType>
auto map(shared_children_t<ValueType, Ratio> children, FunctionType function)
{
    return shared_children_t<ValueType, Ratio>{
        shared_trees(map(*children.ptr, function))
    };
}




/**
 * @brief      Return true if this tree node has a value, false otherwise (in
 *             which case it has Ratio children).
 *
 * @param[in]  tree          The tree to test
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     True if has value, False otherwise.
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
bool has_value(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    return std::holds_alternative<ValueType>(tree.provider);
}




/**
 * @brief      Return the value at this node. Throws an exception if this is not
 *             a leaf node.
 *
 * @param[in]  tree          The tree node
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     The value
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto value(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    return std::get<ValueType>(tree.provider);
}




/**
 * @brief      Return the children of this node. Throws an exception if this is
 *             a leaf node.
 *
 * @param[in]  tree          The tree node
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     The children
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto children(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    if (has_value(tree))
    {
        throw std::out_of_range("bsp::child (tree is a leaf)");
    }
    return std::get<ChildrenType>(tree.provider);
}




/**
 * @brief      Return the i-th child of this (non-leaf) tree node.
 *
 * @param[in]  tree          The tree node
 * @param[in]  i             The child index to get
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     Another tree node
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto child_at(const tree_t<ValueType, ChildrenType, Ratio>& tree, std::size_t i)
{
    if (i >= Ratio)
    {
        throw std::out_of_range("bsp::child_at (index must be <= Ratio)");
    }
    return children(tree)(i);
}




/**
 * @brief      Attach children to a leaf node. The value at this node is
 *             replaced by a set of children, resulting from the attach
 *             function: value -> children.
 *
 * @param[in]  tree             The tree node
 * @param[in]  attach_function  The function: value -> children
 *
 * @tparam     ValueType        The tree value type
 * @tparam     ChildrenType     The tree provider type
 * @tparam     Ratio            The tree ratio
 * @tparam     AttachType       The attach function type
 *
 * @return     Another tree node
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, typename AttachType>
auto attach(tree_t<ValueType, ChildrenType, Ratio> tree, AttachType attach_function)
{
    static_assert(std::is_same_v<std::invoke_result_t<AttachType, ValueType>, ChildrenType>,
        "the attach function must be ValueType -> ChildrenType");

    if (! has_value(tree))
    {
        throw std::invalid_argument("bsp::attach (can only attach leaf nodes)");
    }
    return tree_t<ValueType, ChildrenType, Ratio>{
        attach_function(value(tree))
    };
};




/**
 * @brief      Attach children to all leaf nodes satisfying a predicate. The
 *             value at each descendent leaf node is replaced by a set of child
 *             nodes, resulting from the attach function: value -> children, if
 *             the value of that leaf node satisfies the given predicate.
 *
 * @param[in]  tree             The tree node
 * @param[in]  attach_function  The function: ValueType -> ChildrenType
 * @param[in]  predicate        The predicate: ValueType -> bool
 *
 * @tparam     ValueType        The tree value type
 * @tparam     ChildrenType     The tree provider type
 * @tparam     Ratio            The tree ratio
 * @tparam     AttachType       The attach function type
 * @tparam     PredicateType    The predicate function type
 *
 * @return     Another tree node
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, typename AttachType, typename PredicateType>
auto attach_if(tree_t<ValueType, ChildrenType, Ratio> tree, AttachType attach_function, PredicateType predicate)
{
    static_assert(std::is_same_v<std::invoke_result_t<PredicateType, ValueType>, bool>,
        "the predicate must be ValueType -> bool");

    if (has_value(tree))
    {
        return predicate(value(tree)) ? attach(tree, attach_function) : tree;
    }
    return tree_t<ValueType, ChildrenType, Ratio>{
        map(children(tree), [attach_function, predicate] (auto child)
        {
            return attach_if(child, attach_function, predicate);
        })
    };
}




/**
 * @brief      Replace a leaf node with a node whose children are all leaves.
 *             The children values are determined by the branch function.
 *
 * @param[in]  tree             The tree node
 * @param[in]  branch_function  ValueType -> numeric::array_t<ValueType>
 *
 * @tparam     ValueType        The tree value type
 * @tparam     Ratio            The tree ratio
 * @tparam     BranchType       The type of the branch function
 *
 * @return     Another tree node
 *
 * @note       This function only works on shared trees.
 */
template<typename ValueType, unsigned long Ratio, typename BranchType>
auto branch(shared_tree<ValueType, Ratio> tree, BranchType branch_function)
{
    return attach(tree, [f=branch_function] (auto u) { return shared_values(f(u)); });
};




/**
 * @brief      Branch all leaf nodes satisfying a predicate.
 *
 * @param[in]  tree             The tree node
 * @param[in]  branch_function  ValueType -> numeric::array_t<ValueType>
 * @param[in]  predicate        The predicate: ValueType -> bool
 *
 * @tparam     ValueType        The tree value type
 * @tparam     Ratio            The tree ratio
 * @tparam     BranchType       The type of the branch function
 * @tparam     PredicateType    The predicate function type
 *
 * @return     Another tree node
 *
 * @note       This function only works on shared trees.
 */
template<typename ValueType, unsigned long Ratio, typename BranchType, typename PredicateType>
auto branch_if(shared_tree<ValueType, Ratio> tree, BranchType branch_function, PredicateType predicate)
{
    static_assert(std::is_same_v<std::invoke_result_t<BranchType, ValueType>, numeric::array_t<ValueType, Ratio>>,
        "the branch function must be ValueType -> numeric::array_t<ValueType, Ratio>");

    return attach_if(tree, [f=branch_function] (auto u) { return shared_values(f(u)); }, predicate); 
}




/**
 * @brief      Branch all leaf nodes.
 *
 * @param[in]  tree             The tree node
 * @param[in]  branch_function  ValueType -> numeric::array_t<ValueType>
 *
 * @tparam     ValueType        The tree value type
 * @tparam     Ratio            The tree ratio
 * @tparam     BranchType       The type of the branch function
 *
 * @return     Another tree node
 *
 * @note       This function only works on shared trees.
 */
template<typename ValueType, unsigned long Ratio, typename BranchType>
auto branch_all(shared_tree<ValueType, Ratio> tree, BranchType branch_function)
{
    return branch_if(tree, branch_function, [] (auto) { return true; });
}




/**
 * @brief      Collapse a non-leaf node to a value by recursively mapping a
 *             collapsing function over its descendents.
 *
 * @param[in]  tree              The tree to collapse entirely
 * @param[in]  f                 The collapse function:
 *                               numeric::array_t<ValueType, Ratio> -> ValueType
 *
 * @tparam     ValueType         The tree value type
 * @tparam     ChildrenType      The tree provider type
 * @tparam     Ratio             The tree ratio
 * @tparam     CollapseFunction  The type of the collapse function
 *
 * @return     A value
 *
 * @note       This function differs from reduce in that the function maps an
 *             array of values to a value.
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, typename CollapseFunction>
auto collapse(tree_t<ValueType, ChildrenType, Ratio> tree, CollapseFunction f) -> ValueType
{
    static_assert(std::is_same_v<std::invoke_result_t<CollapseFunction, numeric::array_t<ValueType, Ratio>>, ValueType>,
       "the collapse function must be numeric::array_t<ValueType, Ratio> -> ValueType");

    if (has_value(tree))
    {
        throw std::invalid_argument("bsp::collapse_all (cannot collapse a leaf node)");
    }
    return f(map(numeric::range<Ratio>(), [tree, f] (auto i)
    {
        if (auto child = child_at(tree, i); has_value(child))
        {
            return value(child);
        }
        else
        {
            return collapse(child, f);
        }
    }));
}




/**
 * @brief      Map the leaves of a tree through a function, resulting in a
 *             lazily mapped tree.
 *
 * @param[in]  tree          The tree node
 * @param[in]  function      The function: ValueType -> ResultValueType
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 * @tparam     FunctionType  The mapping function type
 *
 * @return     Another tree node
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, typename FunctionType>
auto map(tree_t<ValueType, ChildrenType, Ratio> tree, FunctionType function)
{
    using result_value_type    = std::invoke_result_t<FunctionType, ValueType>;
    using result_children_type = mapped_children_t<ChildrenType, FunctionType>;
    using result_tree_type     = tree_t<result_value_type, result_children_type, Ratio>;

    if (has_value(tree))
        return result_tree_type{function(value(tree))};

    return result_tree_type{
        result_children_type{children(tree), function}
    };
}




/**
 * @brief      Zip a collection of tree having the same topology to a single
 *             tree of std::tuple.
 *
 * @param[in]  trees         The trees to zip together
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     A tree of tuples.
 */
template<typename... ValueType, typename... ChildrenType, unsigned long Ratio>
auto zip(tree_t<ValueType, ChildrenType, Ratio>... trees)
{
    using result_value_type    = std::tuple<ValueType...>;
    using result_children_type = zipped_children_t<ChildrenType...>;
    using result_tree_type     = tree_t<result_value_type, result_children_type, Ratio>;

    if (all(numeric::array(has_value(trees)...)))
        return result_tree_type{std::tuple(value(trees)...)};

    if (! any(numeric::array(has_value(trees)...)))
        return result_tree_type{result_children_type{std::tuple(children(trees)...)}};

    throw std::invalid_argument("bsp::zip (argument trees have different topology)");
}




/**
 * @brief      Return the number of leaves in a tree
 *
 * @param[in]  tree          The tree whose size is needed
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     The size of the tree
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
std::size_t size(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    return has_value(tree) ? 1 : detail::sum<Ratio>([&] (auto i) { return size(child_at(tree, i)); });
}




/**
 * @brief      Return the tree's maximum depth below this node.
 *
 * @param[in]  tree          The tree whose depth is needed
 * @param[in]  d             Starting depth (mainly for internal use)
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     The tree's depth
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
std::size_t depth(tree_t<ValueType, ChildrenType, Ratio> tree, std::size_t d=0)
{
    if (has_value(tree))
    {
        return d;
    }
    return max(map(numeric::range<Ratio>(), [tree, d] (auto i) { return depth(child_at(tree, i), d + 1); }));
}




/**
 * @brief      Apply a reduction to a tree, such as taking the sum or the
 *             mininum or maximum value.
 *
 * @param[in]  tree          The tree to reduce
 * @param[in]  reducer       The reducer function
 * @param[in]  seed          The seed value
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 * @tparam     Reducer       The reducer function type
 * @tparam     SeedType      The seed type
 *
 * @return     The reduced value
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, typename Reducer, typename SeedType>
SeedType reduce(tree_t<ValueType, ChildrenType, Ratio> tree, Reducer reducer, SeedType seed)
{
    static_assert(std::is_same_v<SeedType, std::invoke_result_t<Reducer, SeedType, ValueType>>,
        "the reducer must be (seed, value) -> seed");

    if (has_value(tree))
    {
        return reducer(seed, value(tree));
    }
    for (std::size_t i = 0; i < Ratio; ++i)
    {
        seed = reduce(child_at(tree, i), reducer, seed);
    }
    return seed;
}




/**
 * @brief      Invoke a side-effect only function on the values of the tree.
 *
 * @param[in]  tree          The tree to invoke the function over
 * @param[in]  function      The function to invoke
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 * @tparam     FunctionType  The function type
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio, typename FunctionType>
void sink(tree_t<ValueType, ChildrenType, Ratio> tree, FunctionType function)
{
    if (has_value(tree))
        function(value(tree));
    else
        for (std::size_t i = 0; i < Ratio; ++i)
            sink(child_at(tree, i), function);
}




/**
 * @brief      Convert a lazy tree (zipped or mapped) into a shared tree.
 *
 * @param[in]  tree          The tree to convert
 *
 * @tparam     ValueType     The tree value type
 * @tparam     ChildrenType  The tree provider type
 * @tparam     Ratio         The tree ratio
 *
 * @return     A shared tree
 */
template<typename ValueType, typename ChildrenType, unsigned long Ratio>
auto to_shared(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    static_assert(! std::is_same_v<ChildrenType, shared_children_t<ValueType, Ratio>>,
        "tree is already shared");

    if (has_value(tree))
    {
        return just<Ratio>(value(tree));
    }
    auto children = numeric::array_t<shared_tree<ValueType, Ratio>, Ratio>();

    for (std::size_t i = 0; i < Ratio; ++i)
    {
        children[i] = to_shared(child_at(tree, i));
    }
    return shared_tree<ValueType, Ratio>{shared_trees(children)};
}




//=============================================================================
template<typename V, typename C, unsigned long R, typename F> auto operator|(tree_t<V, C, R> t, F f) { return f(t); }
template<typename F> auto map  (F f) { return [f] (auto tree) { return map(tree, f); }; }
template<typename F> auto maps (F f) { return [f] (auto tree) { return to_shared(map(tree, f)); }; }
template<typename F> auto mapv (F f) { return [f] (auto tree) { return map(tree, [f] (auto t) { return std::apply(f, t); }); }; }
template<typename F> auto mapvs(F f) { return [f] (auto tree) { return to_shared(map(tree, [f] (auto t) { return std::apply(f, t); })); }; }
inline auto to_shared() { return [] (auto tree) { return to_shared(tree); }; }

} // namespace bsp




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_bsp_tree()
{
    require(size(bsp::just<2>(12)) == 1);
    require(size(bsp::from(1, 2, 3)) == 3);
    require(size(attach(bsp::just<2>(0), [] (int n) { return bsp::shared_values(n + 12, n + 13); })) == 2);
    require(size(branch(bsp::just<2>(0), [] (int n) { return std::array{n + 12, n + 13}; })) == 2);
}

#endif // DO_UNIT_TESTS
