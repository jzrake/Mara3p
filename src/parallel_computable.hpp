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
#include <cstdio>
#include <any>
#include <future>
#include <set>
#include <deque>




/*
 * This file provides a small library to support parallel dependency graph
 * execution, through an idiomatic functional programming model inspired in part
 * by the RxCpp observable interface.
 *
 * It is based around a class template called a `computable', which may be seen
 * as a composable future. Computables can be mapped over and zipped together
 * like any functor. Internally, these operations contruct a directed acyclic
 * graph. Although computables are composed in a semantically functional way,
 * the objects themselves are not immutable. Mapping over a computable gives it
 * a new outgoing edge. Letting a computable go out of scope removes an incoming
 * edge from its upstream computables.
 *
 * A computable implements aspects of std::optional and std::future. When first
 * created, its has_value() method returns false, and calling value() raises an
 * exception. When the computable has no incoming edges but has not yet been
 * evaluated, its eligible() method returns true. When eligible, the computable
 * can be submitted for computation by calling the submit(scheduler) method,
 * where scheduler is a function object () -> std::future<T>. This causes the
 * computable's std::future instance to become valid, and its pending() method
 * to return true. When the computation has finished (the future is ready), the
 * computable's ready() method returns true. Once ready, calling the complete()
 * method invalidates the future (causing pending() to return false), causes
 * value() to return the computation result, and clears the computation
 * std::function object, allowing the release of captured upstream values.
 * Calling complete() also removes all of the computable's outgoing edges,
 * possibly triggering any of its downstream values to become eligible.
 */




//=============================================================================
namespace pr // parallel runtime
{




//=============================================================================
template<typename T>
class unique_deque_t
{
public:

    unique_deque_t() {}
    unique_deque_t(std::initializer_list<T> items) { for (auto item : items) push_back(item); }

    decltype(auto) operator[](std::size_t i) const
    {
        return deque.operator[](i);
    }

    decltype(auto) begin() const
    {
        return deque.begin();
    }

    decltype(auto) end() const
    {
        return deque.end();
    }

    decltype(auto) front() const
    {
        return deque.front();
    }

    decltype(auto) back() const
    {
        return deque.back();
    }

    bool empty() const
    {
        return set.empty();
    }

    auto size() const
    {
        return deque.size();
    }

    auto count(const T& value) const
    {
        return set.count(value);
    }

    void push_back(T&& value)
    {
        if (! count(value))
        {
            deque.push_back(value);
            set.insert(value);
        }
    }

    void push_back(const T& value)
    {
        if (! count(value))
        {
            deque.push_back(value);
            set.insert(value);
        }
    }

    void push_front(T&& value)
    {
        if (! count(value))
        {
            deque.push_front(value);
            set.insert(value);
        }
    }

    void push_front(const T& value)
    {
        if (! count(value))
        {
            deque.push_front(value);
            set.insert(value);
        }
    }

    void pop_front()
    {
        set.erase(deque.front());
        deque.pop_front();
    }

    void pop_back()
    {
        set.erase(deque.back());
        deque.pop_back();
    }

    void erase(const T& value)
    {
        if (count(value))
        {
            set.erase(value);
            deque.erase(std::find(deque.begin(), deque.end(), value));
        }
    }

    void clear()
    {
        deque.clear();
        set.clear();
    }

private:
    std::deque<T> deque;
    std::set<T> set;
};




//=============================================================================
class computable_node_t;
using node_list_t = unique_deque_t<computable_node_t*>;
using async_invoke_t = std::function<std::future<std::any>(std::function<std::any()>)>;




//=============================================================================
inline auto synchronous_execution(std::function<std::any()> computation)
{
    auto task = std::packaged_task<std::any()>(computation);
    task();
    return task.get_future();
}




//=============================================================================
class computable_node_t
{
public:


    //=========================================================================
    computable_node_t(const computable_node_t& other) = delete;
    computable_node_t(const std::any& value) : value(value) {}
    computable_node_t(node_list_t all_incoming)
    {
        for (auto i : all_incoming)
        {
            if (! i->has_value())
            {
                incoming.push_back(i);
            }
        }
        for (auto i : incoming)
        {
            i->outgoing.push_back(this);
        }
    }

    ~computable_node_t()
    {
        for (auto i : incoming)
        {
            i->outgoing.erase(this);
        }
    }

    bool has_value() const
    {
        return value.has_value();
    }

    bool ready() const
    {
        if (pending())
        {
            return future_value.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
        }
        else
        {
            throw std::logic_error("computable_node_t::ready (node is not pending)");
        }
    }

    bool pending() const
    {
        return future_value.valid();
    }

    bool eligible() const
    {
        return incoming.empty() && computation != nullptr;
    }

    void submit(async_invoke_t scheduler)
    {
        if (eligible())
        {
            future_value = scheduler(computation);
        }
        else
        {
            throw std::logic_error("computable_node_t::submit (node is not eligible)");
        }
    }

    void set(std::any new_value)
    {
        if (! pending())
        {
            for (auto o : outgoing)
            {
                o->incoming.erase(this);
            }
            outgoing.clear();
            computation = nullptr;
            value = std::move(new_value);
        }
        else
        {
            throw std::logic_error("computable_node_t::set (node is pending)");
        }
    }

    void complete()
    {
        if (pending())
        {
            set(future_value.get());
        }
        else
        {
            throw std::logic_error("computable_node_t::complete (node is not pending)");
        }
    }

    void compute()
    {
        submit(synchronous_execution);
        complete();
    }

    const auto& get_future() const
    {
        return future_value;
    }

    const auto& outgoing_nodes() const
    {
        return outgoing;
    }

    const auto& incoming_nodes() const
    {
        return incoming;
    }

    auto& name(const char* name)
    {
        node_name = name;
        return *this;
    }

    auto name() const
    {
        return node_name;
    }

    auto& immediate(bool should_be_immediate)
    {
        is_immediate = should_be_immediate;
        return *this;
    }

    auto immediate() const
    {
        return is_immediate;
    }

private:
    //=========================================================================
    std::function<std::any()> computation;
    std::future<std::any> future_value;
    std::any value;
    node_list_t incoming;
    node_list_t outgoing;
    unsigned long node_id;
    bool is_immediate = false;
    const char* node_name = "";
    template<typename T> friend class computable_t;
};




//=============================================================================
template<typename ValueType>
class computable_t
{
public:

    computable_t()
    {
    }

    computable_t(const ValueType& value)
    : g(std::make_shared<computable_node_t>(value))
    {
    }

    computable_t(std::function<ValueType()> computation, node_list_t incoming)
    : g(std::make_shared<computable_node_t>(incoming))
    {
        g->computation = computation;
    }

    const auto& value() const
    {
        if (! g->value.has_value())
        {
            throw std::logic_error("computable_t::value (no value)");
        }
        return std::any_cast<const ValueType&>(g->value);
    }

    bool has_value() const
    {
        return g->has_value();
    }

    bool pending() const
    {
        return g->pending();
    }

    bool eligible() const
    {
        return g->eligible();
    }

    void compute()
    {
        g->compute();
    }

    auto node() const
    {
        return g.get();
    }

    auto& name(const char* name)
    {
        g->name(name);
        return *this;
    }

    auto name() const
    {
        return g->name();
    }

    auto& immediate(bool should_be_immediate)
    {
        g->immediate(should_be_immediate);
        return *this;
    }

    auto immediate() const
    {
        return g->immediate();
    }

private:
    std::shared_ptr<computable_node_t> g;
};

template <typename T>
using computable = computable_t<T>;




//=============================================================================
template<typename ValueType, typename Function>
auto operator|(computable<ValueType> c, Function f)
{
    return f(c);
}

template<typename Function>
auto from(Function f)
{
    using value_type = std::invoke_result_t<Function>;
    return computable<value_type>(f, {});
}

template<typename ValueType>
auto just(ValueType value)
{
    return from([value] () { return value; }).immediate(true);
}

template<typename... ValueType>
auto zip(computable<ValueType>... c)
{
    using value_type = std::tuple<ValueType...>;
    return computable<value_type>([c...] () { return std::tuple(c.value()...); }, {c.node()...}).immediate(true);
}

template<typename ValueType, typename Function>
auto map(computable<ValueType> c, Function f, const char* name=nullptr)
{
    using value_type = std::invoke_result_t<Function, ValueType>;
    return computable<value_type>([c, f] () { return f(c.value()); }, {c.node()}).name(name);
}

template<typename Function>
auto map(Function f, const char* name=nullptr)
{
    return [f, name] (auto c) { return map(c, f, name); };
}

template<typename Function>
auto mapv(Function f, const char* name=nullptr)
{
    return [f, name] (auto c) { return map(c, [f] (auto t) { return std::apply(f, t); }, name); };
}




//=============================================================================
void print_graph(computable_node_t* node, FILE* outfile);




/**
 * @brief      Prints the execution graph for a computable in a format readable
 *             by the graphviz dot utility.
 *
 * @param[in]  c          The computable to print the graph for
 * @param      outfile    The file to write to
 *
 * @tparam     ValueType  The computable value type
 *
 * @note       An example command to generate a PDF from an output file
 *             graph.dot is:
 *
 *             dot -Tpdf -o graph.pdf graph.dot
 */
template<typename ValueType>
void print_graph(pr::computable<ValueType> c, FILE* outfile)
{
    print_graph(c.node(), outfile);
}




//=============================================================================
std::pair<node_list_t, std::deque<unsigned>> topological_sort(computable_node_t* node);




/**
 * @brief      Sort the nodes upstream of and including some computable
 *             according to a stable evaluation order.
 *
 * @param[in]  c          The computable
 *
 * @tparam     ValueType  The computable value type
 *
 * @return     Two sequences of the same length: the first containing the sorted
 *             nodes, and the second containing a sequence of generation of each
 *             node.
 */
template<typename ValueType>
std::pair<node_list_t, std::deque<unsigned>> topological_sort(computable<ValueType> c)
{
    return topological_sort(c.node());
}




//=============================================================================
void compute(computable_node_t* node, async_invoke_t scheduler);




/**
 * @brief      Compute a computable on a (possibly asynchronous) scheduler.
 *
 * @param[in]  c          The computable to compute
 * @param[in]  scheduler  The scheduling function
 *
 * @tparam     ValueType  The computable value type
 */
template<typename ValueType>
void compute(computable<ValueType> c, async_invoke_t scheduler=synchronous_execution)
{
    compute(c.node(), scheduler);
}

} // namespace computable
