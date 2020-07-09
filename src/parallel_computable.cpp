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




#include <map>
#include "parallel_computable.hpp"
#include "parallel_message_queue.hpp"
#include "parallel_thread_pool.hpp"
#include "app_control.hpp"




using namespace mpr;
using node_map_t = std::map<computable_node_t*, unsigned, computable_node_t::comparator_t>;
unsigned long mpr::computable_node_t::last_node_id = 0;




//=============================================================================
template<typename Predicate>
static void collect(computable_node_t* node, Predicate pred, node_set_t& result, node_set_t& passed)
{
    if (pred(node))
    {
        result.insert(node);
    }

    for (auto i : node->incoming_nodes())
    {
        if (! passed.count(i))
        {
            collect(i, pred, result, passed);
        }
    }
    passed.insert(node);
}

template<typename Predicate>
static auto collect(const node_set_t& nodes, Predicate pred)
{
    auto result = node_set_t();
    auto passed = node_set_t();

    for (auto node : nodes)
    {
        collect(node, pred, result, passed);
    }
    return result;
}




//=============================================================================
std::pair<node_set_t, node_set_t> mpr::next_generation(node_set_t eligible, node_set_t completed)
{
    auto next_eligible = node_set_t();
    auto desc_eligible = node_set_t();

    for (auto node : eligible)
    {
        completed.insert(node);

        for (auto o : node->outgoing_nodes())
        {
            desc_eligible.insert(o);
        }
    }

    for (auto node : desc_eligible)
    {
        const auto& i = node->incoming_nodes();

        if (std::all_of(i.begin(), i.end(), [&completed] (auto i) { return completed.count(i); }))
        {
            next_eligible.insert(node);
        }
    }
    return {std::move(next_eligible), std::move(completed)};
}




//=============================================================================
std::vector<node_set_t> mpr::topological_sort(const node_set_t& nodes)
{
    auto current   = collect(nodes, [] (auto n) { return n->incoming_nodes().empty(); });
    auto completed = current;
    auto result    = std::vector<node_set_t>();

    while (! current.empty())
    {
        result.push_back(current);
        std::tie(current, completed) = mpr::next_generation(std::move(current), std::move(completed));
    }
    return result;
}




//=============================================================================
static int closest_upstream_group(computable_node_t* node, const node_map_t& delegations)
{
    auto group_count = std::map<int, unsigned>();

    for (auto i : node->incoming_nodes())
    {
        group_count[delegations.at(i)] += 1;
    }

    auto max_count = 0;
    auto max_group = 0;

    for (auto [group, count] : group_count)
    {
        if (count > max_count)
        {
            max_group = group;
        }
    }
    return max_group;
}




//=============================================================================
static node_map_t delegate_tasks(const std::vector<node_set_t>& generations)
{
    auto delegations = node_map_t();
    auto num_groups = mpi::comm_world().size();

    for (const auto& generation : generations)
    {
        // Treat this generation as having a fixed number of tasks per block.

        if (generation.size() % generations.front().size() == 0)
        {
            auto j = unsigned(0);

            for (auto node : generation)
            {
                delegations.emplace(node, j++ * num_groups / generation.size());
            }
        }
        else // Delegate this generation based on proximity to upstream tasks.
        {
            for (auto node : generation)
            {
                delegations.emplace(node, closest_upstream_group(node, delegations));
            }
        }
    }
    return delegations;
}




//=============================================================================
void mpr::print_graph(std::ostream& stream, const node_set_t& node_set)
{
    auto generations = topological_sort(node_set);
    auto delegations = delegate_tasks(generations);
    stream << "digraph {\n";
    stream << "    ranksep=2;\n";

    for (auto [a, node_group] : delegations)
    {
        for (auto b : a->outgoing_nodes())
        {
            if (delegations.at(a) != delegations.at(b))
            {
                stream << "    " << a->id() << " -> " << b->id() << ";\n";
            }
        }
    }

    for (unsigned group = 0; group < unsigned(mpi::comm_world().size()); ++group)
    {
        stream << "    subgraph cluster_" << group << " {\n";
        stream << "        label = \"Rank " << group << "\";\n";
        stream << "        style = filled;\n";
        stream << "        color = " << group + 1 << ";\n";
        stream << "        colorscheme = spectral9;\n";

        for (auto [a, node_group_a] : delegations)
        {
            for (auto b : a->outgoing_nodes())
            {
                auto node_group_b = delegations.at(b);

                if (node_group_a == group && node_group_b == group)
                {
                     stream << "        " << a->id() << " -> " << b->id() <<  ";\n";
                }
            }
        }
        stream << "    }\n";
    }

    for (auto [a, node_group] : delegations)
    {
        stream
        << "    "
        << a->id()
        << "[shape=" << (a->immediate() ? "ellipse" : "box")
        << ",style=" << (a->immediate() ? "dotted" : "filled")
        << ",label=" << '"' << (std::strlen(a->name()) == 0 ? std::to_string(a->id()) : std::string(a->name())) << '"'
        << "]\n";
    }

    stream << "}\n";
}




//=============================================================================
static execution_monitor_t compute_threaded(const node_set_t& node_set, unsigned num_threads)
{
    auto start       = std::chrono::high_resolution_clock::now();
    auto thread_pool = mara::ThreadPool(num_threads ? num_threads : std::thread::hardware_concurrency());
    auto scheduler   = async_invoke_t(thread_pool.scheduler());
    auto eligible    = collect(node_set, [] (auto n) { return n->eligible(); });
    auto pending     = node_set_t();
    auto completed   = node_set_t();

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
                if (next->eligible())
                {
                    eligible.insert(next);
                }
            }
        }

        completed.clear();
    }

    auto monitor = execution_monitor_t();
    monitor.total_time = (std::chrono::high_resolution_clock::now() - start).count();
    return monitor;
}




//=============================================================================
static execution_monitor_t compute_mpi(const node_set_t& node_set, execution_strategy_t strategy)
{
    using async_load_t = std::pair<computable_node_t*, std::future<std::any>>;
    auto start_time = std::chrono::high_resolution_clock::now();




    // ------------------------------------------------------------------------
    // Collect and sort all upstream nodes, including those already evaluated.
    // Assign the unevaluated tasks to an MPI rank based on divvying the tasks
    // in each generation.
    // ------------------------------------------------------------------------
    auto node_lookup = std::map<unsigned long, computable_node_t*>();

    for (auto [node, group] : delegate_tasks(topological_sort(node_set)))
    {
        if (! node->has_group_number())
        {
            node->assign_to_group(group);
        }
        node_lookup.emplace(node->id(), node);
    }




    auto this_group    = mpi::comm_world().rank();
    auto message_queue = mara::MessageQueue();
    auto thread_pool   = mara::ThreadPool(strategy.num_threads ? strategy.num_threads : std::thread::hardware_concurrency());
    auto eligible      = collect(node_set, [this_group] (auto n) { return n->group() == this_group && n->eligible(); });
    auto delegated     = collect(node_set, [this_group] (auto n) { return n->group() == this_group; });
    auto completed     = collect(node_set, [          ] (auto n) { return n->has_value(); });
    auto pending       = node_set_t();
    auto loading       = std::vector<async_load_t>();


    // ------------------------------------------------------------------------
    // Put each node that is immediately downstream of a node A into the
    // eligible container.
    // ------------------------------------------------------------------------
    auto enqueue_eligible_downstream = [this_group, &eligible] (computable_node_t* node)
    {
        for (auto next : node->outgoing_nodes())
        {
            if (next->group() == this_group && next->eligible())
            {
                eligible.insert(next);
            }
        }
    };

    auto unique_recipients = [this_group] (computable_node_t* node)
    {
        auto result = std::set<int>();

        for (auto next : node->outgoing_nodes())
        {
            if (next->group() != this_group)
            {
                result.insert(next->group());
            }
        }
        return result;
    };

    auto priority = [this_group] (computable_node_t* node)
    {
        auto result = 0;

        for (auto n : node->outgoing_nodes())
        {
            if (n->group() != this_group)
            {
                ++result;
            }
        }
        return result;
    };




    auto iteration = 0;
    auto eval_tick = std::chrono::high_resolution_clock::duration::zero();
    auto eval_dead = std::chrono::high_resolution_clock::duration::zero();




    while (! delegated.empty() || ! message_queue.empty())
    {
        auto start_loop = std::chrono::high_resolution_clock::now();




        // --------------------------------------------------------------------
        // Evaluate each eligible node, if it is delegated to this MPI rank.
        // Enqueue each of the eligible nodes downstream of that node.
        // --------------------------------------------------------------------
        for (auto node : eligible)
        {
            node->submit(thread_pool.scheduler(priority(node)));
            pending.insert(node);
        }




        // --------------------------------------------------------------------
        // Poll the pending tasks. For each one that is ready, call its complete
        // method and place it in the list of completed nodes. Ensure pending
        // tasks have been removed from the list of eligible nodes.
        // --------------------------------------------------------------------
        for (auto node : pending)
        {
            if (node->ready())
            {
                node->complete();
                completed.insert(node);
            }
            eligible.erase(node);
        }




        // --------------------------------------------------------------------
        // Remove completed nodes from the list of pending nodes, and the list
        // of nodes delegated to this MPI rank. Send a future representing the
        // serialized result of each completed node A to each MPI rank that is
        // responsible for any of A's downstream nodes.
        // --------------------------------------------------------------------
        for (auto node : completed)
        {
            pending.erase(node);
            delegated.erase(node);

            if (auto recipients = unique_recipients(node); ! recipients.empty())
            {
                if (strategy.async_serialize)
                {
                    message_queue.push(thread_pool.enqueue(priority(node), [node] () { return node->serialize(); }), recipients, node->id());
                }
                else
                {
                    message_queue.push(node->serialize(), recipients, node->id());
                }
            }
            enqueue_eligible_downstream(node);
        }




        // --------------------------------------------------------------------
        // Receive the serialized result of each node that was completed on
        // another MPI rank, and which has a downstream node delegated to this
        // MPI rank. Save a future deserialize the message, running on a
        // background thread.
        // --------------------------------------------------------------------
        for (auto&& [id, bytes] : message_queue.poll_tags())
        {
            auto node = node_lookup.at(id);

            if (strategy.async_load)
            {
                auto load = thread_pool.enqueue(priority(node), [node] (auto&& b) { return node->get_serializer().deserialize(b); }, std::move(bytes));
                loading.emplace_back(node, std::move(load));
            }
            else
            {
                node->load_from(bytes);
                enqueue_eligible_downstream(node);
            }
        }




        // --------------------------------------------------------------------
        // Get any values that have finished deserializing, set that value to
        // the corresponding nodes, and enqueue any newly eligible nodes.
        // --------------------------------------------------------------------
        if (strategy.async_load)
        {
            for (auto& [node, future_value] : loading)
            {
                if (future_value.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
                {
                    node->set(future_value.get());
                    enqueue_eligible_downstream(node);
                }
            }

            loading.erase(std::remove_if(loading.begin(), loading.end(), [] (const auto& v)
            {
                return ! std::get<1>(v).valid();
            }), loading.end());            
        }




        message_queue.check_outgoing();
        completed.clear();




        auto finish_loop = std::chrono::high_resolution_clock::now();
        eval_tick += (finish_loop - start_loop);
        eval_dead += (finish_loop - start_loop) * pending.empty();
        iteration++;
    }


    // ------------------------------------------------------------------------
    // Assign an empty value to the target nodes to disconnect them from their
    // upstream graph.
    // ------------------------------------------------------------------------
    for (auto node : node_set)
    {
        if (! node->has_value())
        {
            node->set(std::any());
        }
    }
    mpi::comm_world().barrier();

    auto monitor = execution_monitor_t();
    monitor.total_time = (std::chrono::high_resolution_clock::now() - start_time).count();
    monitor.dead_time  = eval_dead.count();
    monitor.eval_time  = eval_tick.count();
    return monitor;
}




//=============================================================================
execution_monitor_t mpr::compute(const node_set_t& node_set, execution_strategy_t strategy)
{
    if (strategy.use_mpi)
    {
        return compute_mpi(node_set, strategy);
    }
    else
    {
        return compute_threaded(node_set, strategy.num_threads);
    }
}
