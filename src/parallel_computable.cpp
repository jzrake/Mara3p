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
        else
        {
            for (auto i : node->incoming_nodes())
            {
                collect(i, pred, result, passed);
            }
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
std::pair<unique_deque_t<computable_node_t*>, std::deque<unsigned>>
mpr::topological_sort(computable_node_t* node)
{
    auto N0 = collect(node, [] (auto n) { return n->has_value() || n->eligible(); });
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
            if (is_or_precedes(o, node))
            {
                const auto& i = o->incoming_nodes();

                if (! std::any_of(i.begin(), i.end(), [&N1] (auto i) { return ! N1.count(i); }))
                {
                    N0.push_back(o);
                    G0.push_back(g + 1);
                }
            }
        }
    }
    return std::pair(N1, G1);
}




//=============================================================================
void mpr::print_graph(computable_node_t* node, FILE* outfile)
{
    auto [nodes, generation] = topological_sort(node);
    auto order = std::map<computable_node_t*, unsigned>();

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
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
        std::fprintf(outfile, "    %ld [shape=box,label=\"%s(%lu %u)\",style=%s];\n",
            i, nodes[i]->name(),
            i, generation[i],
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




#include "parallel_message_queue.hpp"




/**
 * @brief      Compute a node on an MPI communicator.
 *
 * @param      main_node  The node to compute
 *
 * @note       First, assign MPI ranks to the tasks. Some tasks will already be
 *             assigned to an MPI rank. These tasks may or may not have a value.
 *             If they don't have a value, that's because they were delegated to
 *             a different MPI rank. If the task already has a rank assignment,
 *             do not reassign it. Next, begin the main loop, evaluating
 *             eligible nodes as in the single-rank algorithm. When a node is
 *             encountered that has been delegated to a different rank, set its
 *             value to std::any() and move on.
 */
void mpr::compute_mpi(computable_node_t* main_node)
{




    auto this_group    = mpi::comm_world().rank();
    auto message_queue = mara::MessageQueue();
    auto eligible      = collect(main_node, [this_group] (auto n) { return n->group() == this_group && n->eligible(); });
    auto completed     = collect(main_node, [] (auto n) { return n->has_value(); });




    //--------------------------------------------------------------------------
    // Collect and sort all upstream nodes, including those already evaluated.
    // Assign the unevaluated tasks to an MPI rank based on divvying the tasks
    // in each generation.
    // ------------------------------------------------------------------------
    auto [sorted_nodes, generation] = topological_sort(main_node);
    auto delegation = divvy_tasks(sorted_nodes, generation, 4);
    auto order = std::map<computable_node_t*, unsigned>();

    for (std::size_t i = 0; i < sorted_nodes.size(); ++i)
    {
        if (! sorted_nodes[i]->has_group_number())
        {
            sorted_nodes[i]->assign_to_group(delegation[i]);            
        }
        order[sorted_nodes[i]] = i;
    }




    //--------------------------------------------------------------------------
    // Put each node that is immediately downstream of a node A into the
    // eligible container.
    // ------------------------------------------------------------------------
    auto enqueue_eligible_downstream = [main_node, this_group, &eligible] (computable_node_t* node)
    {
        for (auto next : node->outgoing_nodes())
        {
            if (next->group() == this_group && next->eligible() && is_or_precedes(next, main_node))
            {
                eligible.push_back(next);
            }
        }
    };

    auto unique_recipients = [] (computable_node_t* node)
    {
        auto result = std::set<int>();

        for (auto next : node->outgoing_nodes())
        {
            result.insert(next->group());
        }
        return result;
    };




    while (! eligible.empty())
    {




        //----------------------------------------------------------------------
        // Evaluate each eligible node, if it is delegated to this MPI rank.
        // Enqueue each of the nodes downstream of
        //
        // --------------------------------------------------------------------
        for (auto node : eligible)
        {
            node->submit(synchronous_execution);
            node->complete();
            enqueue_eligible_downstream(node);
            completed.push_back(node);
        }




        //---------------------------------------------------------------------
        // Send the serialized result of each completed node A to each MPI rank
        // that is responsible for any of A's downstream nodes.
        // --------------------------------------------------------------------
        for (auto node : completed)
        {
            message_queue.push(node->serialize(), unique_recipients(node), order[node]);
        }




        //----------------------------------------------------------------------
        // Receive the serialized result of each node that was completed on
        // another MPI rank, and which has a downstream node delegated to this
        // MPI rank.
        // --------------------------------------------------------------------
        for (const auto& [index, bytes] : message_queue.poll_tags())
        {
            sorted_nodes[index]->load_from(bytes);
            enqueue_eligible_downstream(sorted_nodes[index]);
        }

        completed.clear();
    }




    //-------------------------------------------------------------------------
    // Assign an empty value to any nodes that do not have one.
    // ------------------------------------------------------------------------
    for (auto node : sorted_nodes)
    {
        if (! node->has_value())
        {
            if (node->group() == this_group)
            {
                throw std::logic_error("mpr::compute_mpi (node delegated to this MPI rank was not evaluated)");
            }
            node->set(std::any());
        }
    }
}
