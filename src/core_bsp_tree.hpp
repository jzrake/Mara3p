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




// bsp := Binary Space Partitioning (kd-tree, quad-tree, oct-tree)
namespace bsp {


using uint = unsigned long;


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
template<typename ValueType, typename ChildrenType, uint Ratio>
struct tree_t
{
    std::variant<ValueType, ChildrenType> provider;
};




//=============================================================================
template<typename ValueType, uint Ratio>
struct shared_children_t
{
    static constexpr uint ratio = Ratio;
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
template<typename ValueType, uint Ratio>
using shared_tree = tree_t<ValueType, shared_children_t<ValueType, Ratio>, Ratio>;

template<uint Ratio, typename ValueType>
auto just(ValueType value)
{
    return shared_tree<ValueType, Ratio>{value};
}

template<typename ValueType, uint Ratio>
auto shared_trees(numeric::array_t<shared_tree<ValueType, Ratio>, Ratio> child_trees)
{
    return shared_children_t<ValueType, Ratio>{
        std::make_shared<decltype(child_trees)>(child_trees)
    };
}

template<typename ValueType, uint Ratio>
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
template<typename ValueType, uint Ratio, typename FunctionType>
auto map(shared_children_t<ValueType, Ratio> children, FunctionType function)
{
    return shared_children_t<ValueType, Ratio>{
        shared_trees(map(*children.ptr, function))
    };
}




/**
 * @brief      Determines if it has value.
 *
 * @param[in]  tree          The tree
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     True if has value, False otherwise.
 */
template<typename ValueType, typename ChildrenType, uint Ratio>
bool has_value(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    return std::holds_alternative<ValueType>(tree.provider);
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree          The tree
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio>
auto value(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    return std::get<ValueType>(tree.provider);
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree          The tree
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio>
auto children(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    if (has_value(tree))
        throw std::out_of_range("bsp::child (tree is a leaf)");
    return std::get<ChildrenType>(tree.provider);
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree          The tree
 * @param[in]  i             { parameter_description }
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio>
auto child_at(const tree_t<ValueType, ChildrenType, Ratio>& tree, std::size_t i)
{
    if (i >= Ratio)
        throw std::out_of_range("bsp::child (index must be <= Ratio)");
    return children(tree)(i);
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree             The tree
 * @param[in]  attach_function  The update function
 *
 * @tparam     ValueType        { description }
 * @tparam     ChildrenType     { description }
 * @tparam     Ratio            { description }
 * @tparam     AttachType       { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio, typename AttachType>
auto attach(tree_t<ValueType, ChildrenType, Ratio> tree, AttachType attach_function)
{
    static_assert(std::is_same_v<std::invoke_result_t<AttachType, ValueType>, ChildrenType>,
        "The attach function must be ValueType -> ChildrenType");

    if (! has_value(tree))
        throw std::invalid_argument("bsp::attach (can only attach leaf nodes)");

    return tree_t<ValueType, ChildrenType, Ratio>{
        attach_function(value(tree))
    };
};




/**
 * @brief      { function_description }
 *
 * @param[in]  tree             The tree
 * @param[in]  attach_function  The attach function
 * @param[in]  predicate        The predicate
 *
 * @tparam     ValueType        { description }
 * @tparam     ChildrenType     { description }
 * @tparam     Ratio            { description }
 * @tparam     AttachType       { description }
 * @tparam     PredicateType    { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio, typename AttachType, typename PredicateType>
auto attach_if(tree_t<ValueType, ChildrenType, Ratio> tree, AttachType attach_function, PredicateType predicate)
{
    static_assert(std::is_same_v<std::invoke_result_t<PredicateType, ValueType>, bool>,
        "The predicate must be ValueType -> bool");

    if (has_value(tree))
        return predicate(value(tree)) ? attach(tree, attach_function) : tree;

    return tree_t<ValueType, ChildrenType, Ratio>{
        map(children(tree), [attach_function, predicate] (auto child)
        {
            return attach_if(child, attach_function, predicate);
        })
    };
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree             The tree
 * @param[in]  branch_function  The branch function
 *
 * @tparam     ValueType        { description }
 * @tparam     Ratio            { description }
 * @tparam     BranchType       { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, uint Ratio, typename BranchType>
auto branch(shared_tree<ValueType, Ratio> tree, BranchType branch_function)
{
    return attach(tree, [f=branch_function] (auto u) { return shared_values(f(u)); });
};




/**
 * @brief      { function_description }
 *
 * @param[in]  tree             The tree
 * @param[in]  branch_function  The branch function
 * @param[in]  predicate        The predicate
 *
 * @tparam     ValueType        { description }
 * @tparam     Ratio            { description }
 * @tparam     BranchType       { description }
 * @tparam     PredicateType    { description }
 * @tparam     <unnamed>        { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, uint Ratio, typename BranchType, typename PredicateType,
typename = std::enable_if_t<std::is_same_v<std::invoke_result_t<BranchType, ValueType>, numeric::array_t<ValueType, Ratio>>>>
auto branch_if(shared_tree<ValueType, Ratio> tree, BranchType branch_function, PredicateType predicate)
{
    return attach_if(tree, [f=branch_function] (auto u) { return shared_values(f(u)); }, predicate); 
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree             The tree
 * @param[in]  branch_function  The branch function
 *
 * @tparam     ValueType        { description }
 * @tparam     Ratio            { description }
 * @tparam     BranchType       { description }
 * @tparam     <unnamed>        { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, uint Ratio, typename BranchType,
typename = std::enable_if_t<std::is_same_v<std::invoke_result_t<BranchType, ValueType>, numeric::array_t<ValueType, Ratio>>>>
auto branch_through(shared_tree<ValueType, Ratio> tree, BranchType branch_function)
{
    return branch_if(tree, branch_function, [] (auto) { return true; });
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree          The tree
 * @param[in]  function      The function
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 * @tparam     FunctionType  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio, typename FunctionType>
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
 * @brief      { function_description }
 *
 * @param[in]  trees         The trees
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename... ValueType, typename... ChildrenType, uint Ratio>
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
 * @param[in]  tree          The tree
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio>
std::size_t size(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    return has_value(tree) ? 1 : detail::sum<Ratio>([&] (auto i) { return size(child_at(tree, i)); });
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree          The tree
 * @param[in]  function      The function
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 * @tparam     FunctionType  { description }
 */
template<typename ValueType, typename ChildrenType, uint Ratio, typename FunctionType>
void sink(tree_t<ValueType, ChildrenType, Ratio> tree, FunctionType function)
{
    if (has_value(tree))
        function(value(tree));
    else
        for (std::size_t i = 0; i < Ratio; ++i)
            sink(child_at(tree, i), function);
}

template<typename ValueType, uint Ratio, typename FunctionType>
void sink(shared_tree<ValueType, Ratio> tree, FunctionType function)
{
    if (has_value(tree))
        function(std::get<ValueType>(tree.provider));
    else
        for (std::size_t i = 0; i < Ratio; ++i)
            sink(std::get<shared_children_t<ValueType, Ratio>>(tree.provider).ptr->operator[](i), function);
}




/**
 * @brief      { function_description }
 *
 * @param[in]  tree          The tree
 *
 * @tparam     ValueType     { description }
 * @tparam     ChildrenType  { description }
 * @tparam     Ratio         { description }
 *
 * @return     { description_of_the_return_value }
 */
template<typename ValueType, typename ChildrenType, uint Ratio>
auto to_shared(tree_t<ValueType, ChildrenType, Ratio> tree)
{
    static_assert(! std::is_same_v<ChildrenType, shared_children_t<ValueType, Ratio>>, "tree is already shared");

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
template<typename V, typename C, uint R, typename F> auto operator|(tree_t<V, C, R> t, F f) { return f(t); }
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
