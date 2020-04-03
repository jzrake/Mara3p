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
#include <any>
#include <future>
#include <ostream>
#include <set>
#include <vector>
#include "app_serial.hpp"




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
 * exception. When a computable's upstream nodes all have values, its eligible()
 * method returns true. When eligible, the computable can be submitted for
 * computation by calling the submit(scheduler) method, where scheduler is a
 * function object () -> std::future<T>. This causes the computable's
 * std::future instance to become valid, and its pending() method to return
 * true. When the computation has finished (the future is ready), the
 * computable's ready() method returns true. Once ready, calling the complete()
 * method invalidates the future (causing pending() to return false), causes
 * value() to return the computation result, and clears the computation
 * std::function object, allowing the release of captured upstream values.
 * Calling complete() also removes all of the computable's outgoing edges,
 * possibly triggering any of its downstream values to become eligible.
 */




//=============================================================================
namespace mpr // massively parallel runtime
{




//=============================================================================
using async_invoke_t = std::function<std::future<std::any>(std::function<std::any()>)>;




//=============================================================================
inline auto synchronous_execution(std::function<std::any()> computation)
{
    auto task = std::packaged_task<std::any()>(computation);
    task();
    return task.get_future();
}




//=============================================================================
struct computable_serializer_t
{
    virtual ~computable_serializer_t() {}
    virtual std::vector<char> serialize(std::any value) const = 0;
    virtual std::any deserialize(const std::vector<char>& bytes) const = 0;
};




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

    using node_set_t = std::set<computable_node_t*, comparator_t>;




    //=========================================================================
    computable_node_t(const computable_node_t& other) = delete;
    computable_node_t(const std::any& any_value) : any_value(any_value), node_id(last_node_id++) {}
    computable_node_t(node_set_t incoming) : incoming(incoming), node_id(last_node_id++)
    {
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

    bool has_value() const
    {
        return any_value.has_value();
    }

    const auto& value() const
    {
        return any_value;
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
        if (has_value())
        {
            return false;
        }
        for (auto node : incoming)
        {
            if (! node->has_value())
            {
                return false;
            }
        }
        return true;
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
            for (auto i : incoming)
            {
                i->outgoing.erase(this);
            }
            incoming.clear();
            computation = nullptr;
            any_value = std::move(new_value);
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

    auto& name(const char* new_name)
    {
        node_name = new_name;
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

    void assign_to_group(int new_group)
    {
        task_group = new_group;
    }

    bool has_group_number()
    {
        return task_group >= 0;
    }

    auto group() const
    {
        return task_group;
    }

    auto serialize() const
    {
        return serializer->serialize(value());
    }

    void load_from(const std::vector<char>& bytes)
    {
        set(serializer->deserialize(bytes));
    }

    const auto& get_serializer() const
    {
        return *serializer;
    }

    auto id() const
    {
        return node_id;
    }

private:
    //=========================================================================
    std::shared_ptr<computable_serializer_t> serializer;
    std::function<std::any()> computation;
    std::future<std::any> future_value;
    std::any any_value;
    node_set_t incoming;
    node_set_t outgoing;
    unsigned long node_id = 0;
    int task_group = -1;
    bool is_immediate = false;
    const char* node_name = "";
    template<typename T> friend class computable_t;
    static unsigned long last_node_id;
};

using node_set_t = computable_node_t::node_set_t;




//=============================================================================
template<typename ValueType>
class computable_t
{
public:

    computable_t()
    {
    }

    computable_t(const ValueType& value)
    {
        g = std::make_shared<computable_node_t>(value);
        g->serializer = std::make_shared<serializer_t>();
    }

    computable_t(std::function<ValueType()> computation, node_set_t incoming)
    {
        g = std::make_shared<computable_node_t>(incoming);
        g->computation = computation;
        g->serializer = std::make_shared<serializer_t>();
    }

    const auto& value() const
    {
        if (! g->has_value())
        {
            throw std::logic_error("computable_t::value (no value)");
        }
        return std::any_cast<const ValueType&>(g->value());
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

    auto& name(const char* new_name)
    {
        g->name(new_name);
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

    //=========================================================================
    struct serializer_t : computable_serializer_t
    {
        std::vector<char> serialize(std::any value) const override
        {
            if constexpr (serial::is_serializable<ValueType>())
            {
                return serial::dumps(std::any_cast<ValueType>(value));
            }
            throw std::runtime_error("mpr::computable (attempt to serialize type lacking serialization traits)");
        }
        std::any deserialize(const std::vector<char> &bytes) const override
        {
            if constexpr (serial::is_serializable<ValueType>())
            {
                return serial::loads<ValueType>(bytes);
            }
            throw std::runtime_error("mpr::computable (attempt to deserialize type lacking serialization traits)");
        }
    };
};

template <typename T>
using computable = computable_t<T>;




//=============================================================================
template<typename... ValueType>
struct zipped_computable_t
{
    std::tuple<computable<ValueType>...> tuple;
};




//=============================================================================
template<typename ValueType, typename Function>
auto operator|(computable<ValueType> c, Function f)
{
    return f(c);
}

template<typename... ValueType, typename Function>
auto operator|(zipped_computable_t<ValueType...> c, Function f)
{
    return f(c);
}

template<typename ValueType>
auto just(ValueType value)
{
    return computable<ValueType>(value).immediate(true);
}

template<typename Function>
auto from(Function f)
{
    using value_type = std::invoke_result_t<Function>;
    return computable<value_type>(f, {});
}

template<typename... ValueType>
auto zip(computable<ValueType>... c)
{
    return zipped_computable_t<ValueType...>{std::tuple(c...)};
}

template<typename ValueType, typename Function>
auto map(computable<ValueType> c, Function f, const char* name="")
{
    using value_type = std::invoke_result_t<Function, ValueType>;
    return computable<value_type>([c, f] () { return std::invoke(f, c.value()); }, {c.node()}).name(name);
}

template<typename... ValueType, typename Function>
auto map(zipped_computable_t<ValueType...> c, Function f, const char* name="")
{
    using value_type = std::invoke_result_t<Function, std::tuple<ValueType...>>;

    return std::apply([name, f] (auto... t)
    {
        return computable<value_type>([f, t...] ()
        {
            return std::invoke(f, std::tuple(t.value()...));
        }, {t.node()...}).name(name);
    }, c.tuple);
}

template<typename Function>
auto map(Function f, const char* name="")
{
    return [f, name] (auto c) { return map(c, f).name(name); };
}

template<typename Function>
auto mapv(Function f, const char* name="")
{
    return [f, name] (auto c) { return map(c, [f] (auto t) { return std::apply(f, t); }).name(name); };
}




/**
 * @brief      Return the generation of nodes that will be eligible following
 *             the completion of the eligible nodes given. This operation is
 *             read-only, it does not change the nodes themselves.
 *
 * @param[in]  eligible   The eligible nodes
 * @param[in]  completed  Nodes already completed
 *
 * @return     A pair of the eligible nodes in the next generation, and the
 *             updated set of completed nodes.
 */
std::pair<node_set_t, node_set_t> next_generation(node_set_t eligible, node_set_t completed);




/**
 * @brief      Sort the graph culminating in the given list of computable nodes
 *             into generations.
 *
 * @param[in]  nodes  The nodes representing the graph to be sorted
 *
 * @return     A vector of the sets, starting with the generation that must be
 *             evaluated first
 *
 * @note       In graphviz, the generations are the rows along which the nodes
 *             are arranged.
 */
std::vector<node_set_t> topological_sort(const node_set_t& nodes);




//=============================================================================
void print_graph(std::ostream& stream, const node_set_t& node_list);




/**
 * @brief      Prints the execution graph for a computable in a format readable
 *             by the graphviz dot utility.
 *
 * @param      outfile    The file to write to
 * @param[in]  cs         The computables to print the graph for
 *
 * @tparam     ValueType  The computable value type
 *
 * @note       An example command to generate a PDF from an output file
 *             graph.dot is:
 *
 *             dot -Tpdf -o graph.pdf graph.dot
 */
template<typename... ValueType>
void print_graph(std::ostream& stream, computable<ValueType>... cs)
{
    print_graph(stream, node_set_t{cs.node()...});
}




template<typename SequenceProtol>
void print_graph_all(std::ostream& stream, SequenceProtol sequence)
{
    auto state = start(sequence);
    auto nodes = node_set_t();

    while (state.has_value())
    {
        nodes.insert(obtain(sequence, state.value()).node());
        state = next(sequence, state.value());
    }
    print_graph(stream, nodes);
}




/**
 * @brief      Struct describing an execution strategy to compute a list of
 *             nodes. MPI execution delegates the work to all ranks
 *             participating in mpi::comm_world. This only guarantees that the
 *             computed value is available on at least one of the MPI ranks. It
 *             is the responsibility of the calling code to maintain ownership
 *             of all the computable_node_t pointers in the node_list argument.
 *             This is done automatically if the corresponding computable
 *             instances exist in some type of container when the compute
 *             function is invoked.
 *
 * @note       The MPI rank to which a node was delegated can be queried by
 *             calling its group() method. Even if a node A was not delegated to
 *             this MPI rank, it may still have a value. This occurs if any
 *             nodes just downstream A are delegated to this MPI rank. In that
 *             case, A will have received the computed value from its delegated
 *             rank.
 */
struct execution_strategy_t
{
    bool use_mpi         = false; // use the MPI execution strategy
    bool async_serialize = false; // serialize data products on a background thread
    bool async_load      = false; // deserialize data products on a background thread
    unsigned num_threads = 1;     // size of the thread pool (use zero for hardware concurrency)
};




/**
 * @brief      Factory functions for common execution strategies
 *
 * @param[in]  num_threads  The number threads
 *
 * @return     The execution strategy.
 */
inline execution_strategy_t multi_threaded_execution(unsigned num_threads)
{
    return {false, false, false, num_threads};
}

inline execution_strategy_t mpi_single_threaded_execution()
{
    return {true, false, false, 1};
}

inline execution_strategy_t mpi_multi_threaded_execution(unsigned num_threads)
{
    return {true, false, false, num_threads};
}




/**
 * @brief      { function_description }
 *
 * @param[in]  node_set  The node set
 * @param[in]  strategy  The strategy
 */
void compute(const node_set_t& node_set, execution_strategy_t strategy);




/**
 * @brief      Calculates the all.
 *
 * @param[in]  sequence        The sequence
 * @param[in]  strategy        The strategy
 *
 * @tparam     SequenceProtol  { description }
 *
 * @return     All.
 */
template<typename SequenceProtol>
auto compute_all(SequenceProtol sequence, execution_strategy_t strategy)
{
    auto state = start(sequence);
    auto nodes = node_set_t();

    while (state.has_value())
    {
        nodes.insert(obtain(sequence, state.value()).node());
        state = next(sequence, state.value());
    }
    compute(nodes, strategy);

    return sequence;
}




/**
 * @brief      Compute a set of computables on using the default hardware
 *             concurrency and no MPI.
 *
 * @param[in]  cs         The computables to compute
 *
 * @tparam     ValueType  The computable value type
 */
template<typename... ValueType>
void compute(computable<ValueType>... cs)
{
    compute({cs.node()...}, multi_threaded_execution(0));
}




/**
 * @brief      { function_description }
 *
 * @param[in]  strategy   The strategy
 * @param[in]  cs         The create struct
 *
 * @tparam     ValueType  { description }
 */
template<typename... ValueType>
void compute(execution_strategy_t strategy, computable<ValueType>... cs)
{
    compute({cs.node()...}, strategy);
}

} // namespace mpr
