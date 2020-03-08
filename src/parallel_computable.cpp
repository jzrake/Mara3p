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
#include "parallel_computable.hpp"




unsigned long pr::computable_node_t::last_node_id = 0;




//=============================================================================
void pr::topological_sort(computable_node_t* node)
{
    auto eligible = std::vector<std::pair<unsigned, computable_node_t*>>();
    auto sorted_v = std::vector<std::pair<unsigned, computable_node_t*>>();
    auto sorted_s = std::set<computable_node_t*>();

    for (auto p : node->primitives())
    {
        eligible.emplace_back(0, p);
    }

    while (! eligible.empty())
    {
        auto [generation, node] = eligible.front();

        eligible.erase(eligible.begin());
        sorted_s.insert(node);
        sorted_v.emplace_back(generation, node);

        for (auto o : node->outgoing_nodes())
        {
            const auto& i = o->incoming_nodes();

            if (! std::any_of(i.begin(), i.end(), [&sorted_s] (auto i) { return ! sorted_s.count(i); }))
            {
                eligible.emplace_back(generation + 1, o);
            }
        }
    }

    for (auto n : sorted_v)
    {
        std::printf("gen %u: %04ld %s\n", n.first, n.second->id(), n.second->name());
    }
}




//=============================================================================
void pr::compute(computable_node_t* main_node, async_invoke_t scheduler)
{
    auto eligible  = main_node->primitives();
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
                if (next->eligible() && next->is_or_precedes(main_node))
                {
                    eligible.insert(next);
                }
            }
        }
        completed.clear();
    }
}
