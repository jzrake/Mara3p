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
unsigned long mpr::computable_node_t::last_node_id = 0;




//=============================================================================
template<typename Predicate>
static void collect(computable_node_t* node, Predicate pred, node_set_t& result, node_set_t& passed)
{
    if (! passed.count(node))
    {
        passed.insert(node);

        if (pred(node))
        {
            result.insert(node);
        }

        for (auto i : node->incoming_nodes())
        {
            collect(i, pred, result, passed);
        }
    }
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
std::pair<std::vector<computable_node_t*>, std::vector<unsigned>> mpr::topological_sort(const node_set_t& nodes)
{
    throw;

    auto N0 = collect(nodes, [] (auto n) { return n->incoming_nodes().empty(); });
    auto N1 = std::vector<computable_node_t*>();

    auto G0 = std::vector<unsigned>(N0.size(), 0);
    auto G1 = std::vector<unsigned>();

    while (! N0.empty())
    {
        auto n = *N0.begin();
        auto g = *G0.begin();

        N0.erase(N0.begin());
        G0.erase(G0.begin());

        N1.push_back(n);
        G1.push_back(g);

        for (auto o : n->outgoing_nodes())
        {
            // const auto& i = o->incoming_nodes();

            if (false) //! std::any_of(i.begin(), i.end(), [&N1] (auto i) { return ! N1.count(i); }))
            {
                N0.insert(o);
                G0.push_back(g + 1);
            }
        }
    }
    return std::pair(N1, G1);
}




//=============================================================================
std::vector<unsigned> mpr::divvy_tasks(const std::vector<computable_node_t*>& nodes, std::vector<unsigned>& generation, unsigned num_groups)
{
    throw;

    auto count_same = [] (const std::vector<unsigned>& items, std::size_t i0)
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

    auto delegation = std::vector<unsigned>();
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
void mpr::print_graph(std::ostream& stream, const node_set_t& node_set)
{
    throw;

    // auto [nodes, generation] = topological_sort(node_set);
    // auto delegation = divvy_tasks(nodes, generation, mpi::comm_world().size());
    // auto order = std::map<computable_node_t*, unsigned>();
    // auto num_groups = (delegation.empty() ? 0 : *std::max_element(delegation.begin(), delegation.end())) + 1;

    // for (std::size_t i = 0; i < nodes.size(); ++i)
    // {
    //     order[nodes[i]] = i;
    // }

    // stream << "digraph {\n";
    // stream << "    ranksep=4;\n";

    // for (auto a : nodes)
    // {
    //     for (auto b : a->outgoing_nodes())
    //     {
    //         if (delegation[order[a]] != delegation[order[b]])
    //         {
    //             stream << "    " << order[a] << " -> " << order[b] << ";\n";
    //         }
    //     }
    // }

    // for (unsigned group = 0; group < num_groups; ++group)
    // {
    //     stream << "    subgraph cluster_" << group << " {\n";
    //     stream << "        label = \"Rank " << group << "\";\n";
    //     stream << "        style = filled;\n";
    //     stream << "        color = " << group + 1 << ";\n";
    //     stream << "        colorscheme = spectral9;\n";

    //     for (auto a : nodes)
    //     {
    //         for (auto b : a->outgoing_nodes())
    //         {
    //             if (delegation[order[a]] == group && delegation[order[b]] == group)
    //             {
    //                  stream << "        " << order[a] << " -> " << order[b] <<  ";\n";
    //             }
    //         }
    //     }
    //     stream << "    }\n";
    // }

    // for (auto a : nodes)
    // {
    //     stream
    //     << "    "
    //     << order[a]
    //     << "[shape=" << (a->immediate() ? "ellipse" : "box")
    //     << ",style=" << (a->immediate() ? "dotted" : "filled")
    //     << ",label=" << '"' << (std::strlen(a->name()) == 0 ? std::to_string(order[a]) : std::string(a->name())) << '"'
    //     << "]\n";
    // }

    // stream << "}\n";
}




//=============================================================================
void mpr::compute(const node_set_t& node_set, unsigned num_threads)
{
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
}




//=============================================================================
void mpr::compute_mpi(const node_set_t& node_set, unsigned num_threads)
{
    throw;

    auto time_start = std::chrono::high_resolution_clock::now();


    // ------------------------------------------------------------------------
    // Collect and sort all upstream nodes, including those already evaluated.
    // Assign the unevaluated tasks to an MPI rank based on divvying the tasks
    // in each generation.
    // ------------------------------------------------------------------------
    auto [sorted_nodes, generation] = topological_sort(node_set);
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

    using async_load_t = std::pair<computable_node_t*, std::future<std::any>>;

    auto this_group    = mpi::comm_world().rank();
    auto message_queue = mara::MessageQueue();
    auto thread_pool   = mara::ThreadPool(num_threads ? num_threads : std::thread::hardware_concurrency());
    auto eligible      = collect(node_set, [this_group] (auto n) { return n->group() == this_group && n->eligible(); });
    auto delegated     = collect(node_set, [this_group] (auto n) { return n->group() == this_group; });
    auto completed     = collect(node_set, [          ] (auto n) { return n->has_value(); });
    auto pending       = std::set<computable_node_t*>();
    auto loading       = std::vector<async_load_t>();


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




    auto async_serialize = true;
    auto async_load      = true;
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
                if (async_serialize)
                {
                    message_queue.push(thread_pool.enqueue(priority(node), [node] () { return node->serialize(); }), recipients, order[node]);
                }
                else
                {
                    message_queue.push(node->serialize(), recipients, order[node]);
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
        for (auto&& [index, bytes] : message_queue.poll_tags())
        {
            auto node = sorted_nodes[index];

            if (async_load)
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
        if (async_load)
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



    auto time_evaluate = std::chrono::high_resolution_clock::now();




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
    mpi::comm_world().invoke([&] () {
        std::printf("\nRank: %d (using %lu threads)\n", mpi::comm_world().rank(), thread_pool.size());
        std::printf("delegate ........... %lf\n", 1e-9 * (time_delegate  - time_start).count());
        std::printf("collect ............ %lf\n", 1e-9 * (time_collect   - time_start).count());
        std::printf("evaluate ........... %lf\n", 1e-9 * (time_evaluate  - time_start).count());
        std::printf("eval iterations .... %d\n", iteration);
        std::printf("eval total time .... %lf\n", 1e-9 * eval_tick.count());
        std::printf("eval dead time ..... %lf (%.1lf%%)\n", 1e-9 * eval_dead.count(), 100.0 * eval_dead / eval_tick);
        std::printf("\n");
    });
}
