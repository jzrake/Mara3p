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




#include <vector>
#include <map>
#include "parallel_computable.hpp"
#include "parallel_message_queue.hpp"
#include "parallel_thread_pool.hpp"
#include "app_control.hpp"




using namespace mpr;




//=============================================================================
template<typename Predicate>
static void collect(computable_node_t* node, Predicate pred, node_list_t& result, node_list_t& passed)
{
    if (! passed.count(node))
    {
        passed.push_back(node);

        if (pred(node))
        {
            result.push_back(node);
        }

        for (auto i : node->incoming_nodes())
        {
            collect(i, pred, result, passed);
        }
    }
}

template<typename Predicate>
static auto collect(computable_node_t* node, Predicate pred)
{
    auto result = node_list_t();
    auto passed = node_list_t();
    collect(node, pred, result, passed);
    return result;
}

template<typename Predicate>
static auto collect(const node_list_t& nodes, Predicate pred)
{
    auto result = node_list_t();
    auto passed = node_list_t();

    for (auto node : nodes)
    {
        collect(node, pred, result, passed);
    }
    return result;
}

static bool is_or_precedes(computable_node_t* node, computable_node_t* other)
{
    if (other == node || node->outgoing_nodes().count(other))
    {
        return true;
    }
    for (auto o : node->outgoing_nodes())
    {
        if (is_or_precedes(o, other))
        {
            return true;
        }
    }
    return false;
}




//=============================================================================
std::pair<node_list_t, std::deque<unsigned>> mpr::topological_sort(const node_list_t& nodes)
{
    auto N0 = collect(nodes, [] (auto n) { return n->incoming_nodes().empty(); });
    auto N1 = unique_deque_t<computable_node_t*>();

    auto G0 = std::deque<unsigned>(N0.size(), 0);
    auto G1 = std::deque<unsigned>();

    while (! N0.empty())
    {
        auto n = N0.front();
        auto g = G0.front();

        N0.pop_front();
        G0.pop_front();

        N1.push_back(n);
        G1.push_back(g);

        for (auto o : n->outgoing_nodes())
        {
            const auto& i = o->incoming_nodes();

            if (! std::any_of(i.begin(), i.end(), [&N1] (auto i) { return ! N1.count(i); }))
            {
                N0.push_back(o);
                G0.push_back(g + 1);
            }
        }
    }
    return std::pair(N1, G1);
}




//=============================================================================
std::deque<unsigned> mpr::divvy_tasks(const node_list_t& nodes, std::deque<unsigned>& generation, unsigned num_groups)
{
    auto count_same = [] (const std::deque<unsigned>& items, std::size_t i0)
    {
        for (std::size_t i = i0; i < items.size(); ++i)
        {
            if (items[i] != items[i0])
            {
                return i - i0;
            }
        }
        return items.size() - i0;
    };

    auto delegation = std::deque<unsigned>();
    auto gen_start = std::size_t(0);

    while (gen_start < nodes.size())
    {
        auto gen_size = count_same(generation, gen_start);

        for (std::size_t r = 0; r < num_groups; ++r)
        {
            auto start = gen_start + (r + 0) * gen_size / num_groups;
            auto final = gen_start + (r + 1) * gen_size / num_groups;

            for (std::size_t i = start; i < final; ++i)
            {
                delegation.push_back(r);
            }
        }
        gen_start += gen_size;
    }
    return delegation;
}




//=============================================================================
void mpr::print_graph(FILE* outfile, const node_list_t& node_list)
{
    auto [nodes, generation] = topological_sort(node_list);
    auto delegation = divvy_tasks(nodes, generation, mpi::comm_world().size());
    auto order = std::map<computable_node_t*, unsigned>();

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        if (! nodes[i]->has_group_number())
        {
            nodes[i]->assign_to_group(delegation[i]);            
        }
        order[nodes[i]] = i;
    }

    std::fprintf(outfile, "digraph {\n");
    std::fprintf(outfile, "    ranksep=1.0;\n");

    for (auto n : nodes)
    {
        for (auto o : n->outgoing_nodes())
        {
            std::fprintf(outfile, "    %u -> %u;\n", order[n], order[o]);
        }
    }

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        std::fprintf(outfile, "    %ld [shape=box,label=\"%s(%lu: %d)\",style=%s];\n",
            i, nodes[i]->name(),
            i, nodes[i]->group(),// generation[i],
            nodes[i]->immediate() ? "dotted" : "filled");
    }
    std::fprintf(outfile, "}\n");
}




//=============================================================================
void mpr::compute(computable_node_t* main_node, async_invoke_t scheduler)
{
    auto eligible  = collect(main_node, [] (auto n) { return n->eligible(); });
    auto pending   = node_list_t();
    auto completed = node_list_t();

    while (! eligible.empty() || ! pending.empty())
    {
        for (auto node : eligible)
        {
            node->submit(node->immediate() ? synchronous_execution : scheduler);
            pending.push_back(node);
        }

        for (auto node : pending)
        {
            if (node->ready())
            {
                completed.push_back(node);
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
                if (next->eligible() && is_or_precedes(next, main_node))
                {
                    eligible.push_back(next);
                }
            }
        }
        completed.clear();
    }
}




//=============================================================================
void mpr::compute_mpi(const node_list_t& node_list, unsigned num_threads)
{
    auto time_start = std::chrono::high_resolution_clock::now();


    // ------------------------------------------------------------------------
    // Collect and sort all upstream nodes, including those already evaluated.
    // Assign the unevaluated tasks to an MPI rank based on divvying the tasks
    // in each generation.
    // ------------------------------------------------------------------------
    auto [sorted_nodes, generation] = topological_sort(node_list);
    auto delegation = divvy_tasks(sorted_nodes, generation, mpi::comm_world().size());
    auto order = std::map<computable_node_t*, unsigned>();

    for (std::size_t i = 0; i < sorted_nodes.size(); ++i)
    {
        if (! sorted_nodes[i]->has_group_number())
        {
            sorted_nodes[i]->assign_to_group(delegation[i]);            
        }
        order[sorted_nodes[i]] = i;
    }


    auto time_delegate = std::chrono::high_resolution_clock::now();



    auto this_group    = mpi::comm_world().rank();
    auto message_queue = mara::MessageQueue();
    auto thread_pool   = mara::ThreadPool(num_threads ? num_threads : std::thread::hardware_concurrency());
    auto eligible      = collect(node_list, [this_group] (auto n) { return n->group() == this_group && n->eligible(); });
    auto delegated     = collect(node_list, [this_group] (auto n) { return n->group() == this_group; }).item_set();
    auto completed     = collect(node_list, [          ] (auto n) { return n->has_value(); }).item_set();
    auto pending       = std::set<computable_node_t*>();



    auto time_collect = std::chrono::high_resolution_clock::now();



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
                eligible.push_back(next);
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


    auto iteration = 0;
    auto eval_tick = std::chrono::high_resolution_clock::duration::zero();
    auto eval_dead = std::chrono::high_resolution_clock::duration::zero();
    auto eval_dump = std::chrono::high_resolution_clock::duration::zero();
    auto eval_load = std::chrono::high_resolution_clock::duration::zero();


    while (! delegated.empty())
    {
        auto start_loop = std::chrono::high_resolution_clock::now();



        // --------------------------------------------------------------------
        // Evaluate each eligible node, if it is delegated to this MPI rank.
        // Enqueue each of the eligible nodes downstream of that node.
        // --------------------------------------------------------------------
        for (auto node : eligible)
        {
            node->submit(thread_pool.scheduler());
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
        // of nodes delegated to this MPI rank. Send the serialized result of
        // each completed node A to each MPI rank that is responsible for any of
        // A's downstream nodes.
        // --------------------------------------------------------------------
        for (auto node : completed)
        {
            pending.erase(node);
            delegated.erase(node);

            if (auto recipients = unique_recipients(node); ! recipients.empty())
            {
                auto [bytes, ticks] = control::invoke_timed(std::mem_fn(&computable_node_t::serialize), *node);
                eval_dump += ticks;

                // auto bytes = node->serialize();
                message_queue.push(bytes, recipients, order[node]);
            }
            enqueue_eligible_downstream(node);
        }




        // --------------------------------------------------------------------
        // Receive the serialized result of each node that was completed on
        // another MPI rank, and which has a downstream node delegated to this
        // MPI rank.
        // --------------------------------------------------------------------
        for (const auto& [index, bytes] : message_queue.poll_tags())
        {
            auto node = sorted_nodes[index];
            eval_load += control::invoke_timed(std::mem_fn(&computable_node_t::load_from), *node, bytes);
            // node->load_from(bytes);
            enqueue_eligible_downstream(node);
        }

        completed.clear();


        auto finish_loop = std::chrono::high_resolution_clock::now();


        eval_tick += (finish_loop - start_loop);
        eval_dead += (finish_loop - start_loop) * pending.empty();
        iteration++;
    }



    auto time_evaluate = std::chrono::high_resolution_clock::now();



    // ------------------------------------------------------------------------
    // Assign an empty value to the target nodes to disconnect them from their
    // upstream graph.
    // ------------------------------------------------------------------------
    for (auto node : node_list)
    {
        if (! node->has_value())
        {
            node->set(std::any());
        }
    }



    mpi::comm_world().barrier();



    mpi::comm_world().invoke([&] () {
        std::printf("\nRank: %d (using %lu threads)\n", mpi::comm_world().rank(), thread_pool.size());
        std::printf("delegate ........... %lf\n", 1e-9 * (time_delegate  - time_start).count());
        std::printf("collect ............ %lf\n", 1e-9 * (time_collect   - time_start).count());
        std::printf("evaluate ........... %lf\n", 1e-9 * (time_evaluate  - time_start).count());
        std::printf("eval iterations .... %d\n", iteration);
        std::printf("eval total time .... %lf\n", 1e-9 * eval_tick.count());
        std::printf("eval dead time ..... %lf (%.1lf%%)\n", 1e-9 * eval_dead.count(), 100.0 * eval_dead / eval_tick);
        std::printf("eval dump time ..... %lf (%.1lf%%)\n", 1e-9 * eval_tick.count(), 100.0 * eval_dump / eval_tick);
        std::printf("eval load time ..... %lf (%.1lf%%)\n", 1e-9 * eval_dead.count(), 100.0 * eval_load / eval_tick);
        std::printf("\n");
    });
}
