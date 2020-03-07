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




using async_invoke_t = std::function<std::future<std::any>(std::function<std::any()>)>;




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
    struct comparator_t
    {
        bool operator()(computable_node_t* a, computable_node_t* b) const
        {
            return a->node_id < b->node_id;
        }
    };

    using set_t = std::set<computable_node_t*, comparator_t>;


    //=========================================================================
    computable_node_t(const computable_node_t& other) = delete;
    computable_node_t(const std::any& value) : value(value), node_id(++last_node_id)
    {
    }

    computable_node_t(set_t all_incoming) : node_id(++last_node_id)
    {
        for (auto i : all_incoming)
        {
            if (! i->has_value())
            {
                incoming.insert(i);
            }
        }
        for (auto i : incoming)
        {
            i->outgoing.insert(this);
        }
    }

    ~computable_node_t()
    {
        for (auto i : incoming)
        {
            i->outgoing.erase(this);
        }
    }

    auto id() const
    {
        return node_id;
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

    const auto& first_incoming() const
    {
        return *incoming.begin();
    }

    void primitives(set_t& result)
    {
        if (! result.count(this))
        {
            if (eligible())
            {
                result.insert(this);
            }
            else
            {
                for (auto i : incoming)
                {
                    i->primitives(result);
                }
            }
        }
    }

    auto primitives()
    {
        auto result = set_t();
        primitives(result);
        return result;
    }

    bool is_or_precedes(computable_node_t* other) const
    {
        if (other == this || outgoing.count(other))
        {
            return true;
        }
        for (auto o : outgoing)
        {
            if (o->is_or_precedes(other))
            {
                return true;
            }
        }
        return false;
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
    set_t incoming;
    set_t outgoing;
    unsigned long node_id;
    bool is_immediate = false;
    const char* node_name = "";
    static unsigned long last_node_id;
    template<typename T> friend class computable_t;
};

unsigned long computable_node_t::last_node_id = 0;




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

    computable_t(std::function<ValueType()> computation, computable_node_t::set_t incoming)
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
auto operator|(computable_t<ValueType> c, Function f)
{
    return f(c);
}

template<typename Function>
auto from(Function f)
{
    using value_type = std::invoke_result_t<Function>;
    return computable_t<value_type>(f, {});
}

template<typename ValueType>
auto just(ValueType value)
{
    return from([value] () { return value; }).immediate(true);
}

template<typename... ValueType>
auto zip(computable_t<ValueType>... c)
{
    using value_type = std::tuple<ValueType...>;
    return computable_t<value_type>([c...] () { return std::tuple(c.value()...); }, {c.node()...}).immediate(true);
}

template<typename ValueType, typename Function>
auto map(computable_t<ValueType> c, Function f, const char* name=nullptr)
{
    using value_type = std::invoke_result_t<Function, ValueType>;
    return computable_t<value_type>([c, f] () { return f(c.value()); }, {c.node()}).name(name);
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
void print_recurse(FILE* outfile, pr::computable_node_t* node, pr::computable_node_t::set_t& already_printed)
{
    if (! already_printed.count(node))
    {
        already_printed.insert(node);

        for (auto o : node->outgoing_nodes())
        {
            std::fprintf(outfile, "    %ld -> %ld;\n", node->id(), o->id());
            print_recurse(outfile, o, already_printed);
        }
    }
};




/**
 * @brief      Prints the execution graph for a computable in a format readable
 *             by the graphviz dot utility.
 *
 * @param      outfile  The file to write to
 * @param[in]  c        The computable to print the graph for
 *
 * @tparam     T        The computable value type
 *
 * @note       An example command to generate a PDF from an output file
 *             graph.dot is:
 *             
 *             dot -Tpdf -o graph.pdf graph.dot
 */
template<typename T>
void print_graph(FILE* outfile, pr::computable_t<T> c)
{
    auto already_printed = pr::computable_node_t::set_t();
    std::fprintf(outfile, "digraph {\n");
    std::fprintf(outfile, "    ranksep=1.0;\n");

    for (auto node : c.node()->primitives())
    {
        print_recurse(outfile, node, already_printed);        
    }
    for (auto node : already_printed)
    {
        std::fprintf(outfile, "    %ld [shape=box,label=\"%s\",style=%s];\n",
            node->id(),
            node->name(),
            node->immediate() ? "dotted" : "filled");
    }
    std::fprintf(outfile, "}\n");
}




/**
 * @brief      Compute a computable on a (possibly asynchronous) scheduler.
 *
 * @param[in]  computable  The computable to compute
 * @param[in]  scheduler   The scheduling function
 *
 * @tparam     ValueType   The computable value type
 */
template<typename ValueType>
void compute(computable_t<ValueType> computable, async_invoke_t scheduler=synchronous_execution)
{
    auto eligible  = computable.node()->primitives();
    auto pending   = computable_node_t::set_t();
    auto completed = computable_node_t::set_t();

    while (! eligible.empty() || ! pending.empty())
    {
        for (auto node : eligible)
        {
            node->submit(node->immediate() ? synchronous_execution : scheduler);
            pending.insert(node);
        }

        for (auto node : pending)
        {
            if (node->ready())
            {
                completed.insert(node);
            }
            eligible.erase(node);
        }

        for (auto node : completed)
        {
            auto outgoing = node->outgoing_nodes();

            pending.erase(node);
            node->complete();

            for (auto next : outgoing)
            {
                if (next->eligible() && next->is_or_precedes(computable.node()))
                {
                    eligible.insert(next);
                }
            }
        }
        completed.clear();
    }
}

} // namespace computable
