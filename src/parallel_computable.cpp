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




//=============================================================================
static void primitives(pr::computable_node_t* node, pr::node_list_t& result, pr::node_list_t& passed)
{
    if (! passed.count(node))
    {
        passed.push_back(node);

        if (node->eligible())
        {
            result.push_back(node);
        }
        else
        {
            for (auto i : node->incoming_nodes())
            {
                primitives(i, result, passed);
            }
        }
    }
}

static auto primitives(pr::computable_node_t* node)
{
    auto result = pr::node_list_t();
    auto passed = pr::node_list_t();
    primitives(node, result, passed);
    return result;
}




//=============================================================================
std::pair<pr::unique_deque_t<pr::computable_node_t*>, std::deque<unsigned>>
pr::topological_sort(computable_node_t* node)
{
    auto N0 = primitives(node);
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
            if (o->is_or_precedes(node))
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
void pr::print_graph(computable_node_t* node, FILE* outfile)
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
void pr::compute(computable_node_t* main_node, async_invoke_t scheduler)
{
    auto eligible  = primitives(main_node);
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
                if (next->eligible() && next->is_or_precedes(main_node))
                {
                    eligible.push_back(next);
                }
            }
        }
        completed.clear();
    }
}
