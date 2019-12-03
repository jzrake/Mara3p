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
#include <algorithm>
#include <future>
#include <map>
#include <vector>




//=============================================================================
namespace mara {
    template<typename, typename>
    class DependencyGraph;
}




//=============================================================================
template<typename KeyType, typename ValueType>
class mara::DependencyGraph
{
public:




    //=========================================================================
    using key_type             = KeyType;
    using value_type           = ValueType;
    using mapping_type         = std::function<value_type(std::vector<value_type>)>;
    using key_vector_type      = std::vector<key_type>;
    using rule_type            = std::pair<mapping_type, key_vector_type>;




    /**
     * @brief      Return true if the addition of the given rule would create a
     *             dependency cycle in the graph. This checks for whether any of
     *             the expression's symbols are downstream of its key.
     *
     * @param[in]  key            The key identifying the rule
     * @param[in]  argument_keys  The names of the products required for this
     *                            rule to be evaluated
     *
     * @tparam     KeyTypes       The types of the upstream keys (must all be
     *                            convertable to key_type)
     *
     * @return     True or false
     */
    template<typename... KeyTypes>
    bool cyclic(key_type key, KeyTypes... argument_keys)
    {
        auto dependents = referencing(key);

        for (const auto& s : key_vector_type{argument_keys...})
        {
            if (dependents.count(s) || s == key)
            {
                return true;
            }
        }
        return false;
    }




    /**
     * @brief      Insert an evaluation rule into this graph. The entry must not
     *             replace a rule that already exists.
     *
     * @param[in]  key            The key identifying the rule
     * @param[in]  mapping        The function evaluating the rule
     * @param[in]  argument_keys  The names of the products required for this
     *                            rule to be evaluated
     *
     * @tparam     KeyTypes       The types of the upstream keys (must all be
     *                            convertable to key_type)
     */
    template<typename... KeyTypes>
    void insert_rule(key_type key, mapping_type mapping, KeyTypes... argument_keys)
    {
        if (contains_rule(key))
            throw std::invalid_argument("DependencyGraph::insert_rule (rule already exists)");

        if (cyclic(key, argument_keys...))
            throw std::invalid_argument("DependencyGraph::insert_rule (rule would create a dependency cycle)");

        rules[key] = std::pair(mapping, key_vector_type{argument_keys...});
        downstream[key] = downstream_keys(key);

        for (const auto& [k, r] : rules)
        {
            if (contains(upstream_keys(key), k))
            {
                downstream[k].insert(key);
            }
        }
    }




    /**
     * @brief      Insert a data product into the graph under a given key. The
     *             key must already exist as a rule, but the data product must
     *             not be already evaluated, or pending evaluation. This method
     *             is most likely to be used to insert a data product that was
     *             fetched from a graph on a remote process.
     *
     * @param[in]  item  A std::pair of the key and data product value
     */
    void insert_product(std::pair<key_type, value_type> item)
    {
        if (! contains_rule(item.first))
            throw std::invalid_argument("DependencyGraph::insert_product (rule not defined)");

        if (is_pending(item.first))
            throw std::invalid_argument("DependencyGraph::insert_product (evaluation is pending)");

        if (contains_product(item.first))
            throw std::invalid_argument(std::string("DependencyGraph::insert_product (already evaluated) ") + item.first + " on " + std::to_string(mpi::comm_world().rank()));

        products.insert(item);
    }




    /**
     * @brief      Trigger the asynchronous evaluation of the rule under the
     *             given key, on a particular scheduler. The key must already
     *             exist as a rule, but the data product must not be already
     *             evaluated, or pending evaluation. The scheduler is probably a
     *             thread pool, but the only requirement is that it provides an
     *             enqueue method that returns a std::future<value_type>.
     *
     * @param[in]  key        The key if the data product to evaluate
     * @param      scheduler  The scheduler to use
     *
     * @tparam     Scheduler  The type of the scheduler
     */
    template<typename Scheduler>
    void evaluate_rule(key_type key, Scheduler& scheduler)
    {
        if (! rules.count(key) || products.count(key) || pending_products.count(key))
            throw std::invalid_argument("DependencyGraph::evaluate_rule");

        pending_products[key] = scheduler.enqueue(rules.at(key).first, argument_values(key));
    }




    /**
     * @brief      Check on the status of data products whose evaluation is
     *             currently pending, and move them to the internal map of
     *             completed items.
     *
     * @param[in]  timeout  The amount of time to wait for pending items to
     *                      finish
     *
     * @return     A map of the newly completed items
     */
    std::map<key_type, value_type> poll(std::chrono::duration<long> timeout=std::chrono::seconds(0))
    {
        auto completed = std::map<key_type, value_type>();

        for (auto& [key, future] : pending_products)
        {
            if (future.wait_for(timeout) == std::future_status::ready)
            {
                products[key] = completed[key] = future.get();
            }
        }

        for (const auto& item : completed)
        {
            pending_products.erase(item.first);
        }
        return completed;
    }




    /**
     * @brief      Determine whether the graph has a rule defined under the
     *             given key.
     *
     * @param[in]  key   The key to check for
     *
     * @return     True or false
     */
    bool contains_rule(key_type key) const
    {
        return rules.find(key) != rules.end();
    }




    /**
     * @brief      Determine whether the graph has a product evaluated under the
     *             given key.
     *
     * @param[in]  key   The key to check for
     *
     * @return     True or false
     */
    bool contains_product(key_type key) const
    {
        return products.find(key) != products.end();
    }




    /**
     * @brief      Determine whether an evaluation of the given rule is
     *             currently pending.
     *
     * @param[in]  key   The key to check for
     *
     * @return     True or false
     */
    bool is_pending(key_type key) const
    {
        return pending_products.find(key) != pending_products.end();
    }




    /**
     * @brief      Determines if a given rule is eligible to be evaluated,
     *             meaning that its data product is neither pending nor
     *             complete, all of its argument (upstream) keys are complete,
     *             and that the graph on this process is responsible for
     *             evaluating it. The supplied predicate must return false if
     *             the rule has been delegated to a remote process.
     *
     * @param[in]  key                 The key to check for
     * @param[in]  is_responsible_for  The predicate
     *
     * @return     True if this rule is eligible to be evaluated
     */
    bool is_eligible(key_type key, std::function<bool(key_type)> is_responsible_for=true_predicate()) const
    {
        if (is_responsible_for(key) && ! contains_product(key) && ! is_pending(key))
        {
            for (const auto& argument_key : upstream_keys(key))
            {
                if (! contains_product(argument_key))
                {
                    return false;
                }
            }
            return true;
        }
        return false;
    }




    /**
     * @brief      Return all the rules that are currently eligible for
     *             evaluation.
     *
     * @param[in]  is_responsible_for  The predicate, which must returm false if
     *                                 the given rulehas been delegated to a
     *                                 remote process.
     *
     * @return     A std::vector of keys
     */
    std::vector<key_type> eligible_rules(std::function<bool(key_type)> is_responsible_for=true_predicate()) const
    {
        auto result = std::vector<key_type>();
        
        for (const auto& [key, rule] : rules)
        {
            if (is_eligible(key, is_responsible_for))
            {
                result.push_back(key);
            }
        }
        return result;
    }




    /**
     * @brief      Count the number of unevaluated rules in the graph which
     *             satisfy the given predicate.
     *
     * @param[in]  is_responsible_for  The predicate (must return false if the
     *                                 rule under the given key is delegated to
     *                                 a remote process).
     *
     * @return     The number of unevaluated rules
     */
    unsigned count_unevaluated(std::function<bool(key_type)> is_responsible_for=true_predicate()) const
    {
        unsigned count = 0;

        for (const auto& [key, rule] : rules)
        {
            if (! contains_product(key) && is_responsible_for(key))
            {
                ++count;
            }
        }
        return count;
    }




    /**
     * @brief      Return the immediate dependencies of the rule defined under a
     *             key.
     *
     * @param[in]  key   The key to check for
     *
     * @return     A std::vector of keys
     */
    const std::vector<key_type>& upstream_keys(key_type key) const
    {
        return rules.at(key).second;
    }




    /**
     * @brief      Return the keys of rules immediately downstream of this one
     *             (those for which the given key is immediately upstream).
     *
     * @param[in]  key   The key to check for
     *
     * @return     A std::vector, containing the keys that are immediately
     *             downstream of the one given
     */
    std::set<key_type> downstream_keys(key_type key) const
    {
        if (downstream.count(key))
        {
            return downstream.at(key);            
        }

        auto result = std::set<key_type>();

        for (const auto& [k, r] : rules)
        {
            if (contains(upstream_keys(k), key))
            {
                result.insert(k);
            }
        }
        return result;
    }




    /**
     * @brief      Return the keys that reference (directly or indirectly) the
     *             given key.
     *
     * @param[in]  key   The key to check for
     *
     * @return     A std::set of keys
     */
    std::set<key_type> referencing(key_type key) const
    {
        auto result = downstream_keys(key);

        for (const auto& k : downstream_keys(key))
        {
            for (const auto& m : referencing(k))
            {
                result.insert(m);
            }
        }
        return result;
    }




    /**
     * @brief      Return the data products associated with the rules upstream
     *             of the one given. This function throws if the rule defined
     *             under the given key has unevaluated upstream rules.
     *
     * @param[in]  key   The key whose argument values are required
     *
     * @return     The data products of the arguments for the given key
     */
    std::vector<value_type> argument_values(key_type key) const
    {
        auto result = std::vector<value_type>();

        for (const auto& argument_key : upstream_keys(key))
        {
            result.push_back(products.at(argument_key));
        }
        return result;
    }




    /**
     * @brief      Return the data product at the given key. This function
     *             throws if the rule has not yet been evaluated.
     *
     * @param[in]  key   The key whose product is needed
     *
     * @return     The data product
     */
    const value_type& product_at(key_type key) const
    {
        return products.at(key);
    }




    /**
     * @brief      Return a vector of all the keys in the graph.
     *
     * @return     A std::vector<key_type>
     */
    std::vector<key_type> keys() const
    {
        auto result = std::vector<key_type>();

        for (const auto& item : rules)
        {
            result.push_back(item.first);
        }
        return result;
    }




    /**
     * @brief      Return a status code indicating whether a rule under the
     *             given key is non-existent (0), defined but neither evaluated
     *             nor pending (1), pending evaluation (2), or completed (3).
     *
     * @param[in]  key   The key to check for
     *
     * @return     A status code
     */
    int status(key_type key) const
    {
        if (contains_product(key)) return 3;
        if (is_pending      (key)) return 2;
        if (contains_rule   (key)) return 1;
        return 0;
    }




private:
    //=========================================================================
    static std::function<bool(key_type)> true_predicate()
    {
        return [] (auto) { return true; };
    }

    static bool contains(const std::vector<key_type>& v, const key_type& k)
    {
        return std::find(v.begin(), v.end(), k) != v.end();
    }


    //=========================================================================
    std::map<key_type, rule_type>               rules;
    std::map<key_type, value_type>              products;
    std::map<key_type, std::future<value_type>> pending_products;
    std::map<key_type, std::set<key_type>>      downstream;
};
