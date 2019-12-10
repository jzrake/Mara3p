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




#include <iomanip>
#include <sstream>
#include <variant>
#include "app_hdf5.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_geometric.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_rational.hpp"
#include "app_ui.hpp"
#include "parallel_dependency_graph.hpp"
#include "parallel_thread_pool.hpp"
#include "scheme_mhd_rules.hpp"




using namespace mhd_scheme_v2;
using namespace mhd_rules;
using DependencyGraph = mara::DependencyGraph<product_identifier_t, product_t>;




//=============================================================================
namespace h5 {

void write(const Group& group, std::string name, const product_t& product)
{
    std::visit([&group, name] (auto p) { write(group, name, p); }, product);
}

std::string legalize(std::string s)
{
    return seq::view(s)
    | seq::map([] (auto c) { return c == '/' ? '@' : c; })
    | seq::remove_if([] (auto c) { return c == ' '; })
    | seq::to<std::basic_string>();
}

}




//=============================================================================
void represent_graph_item(std::ostream& os, const DependencyGraph& graph, product_identifier_t key)
{
    auto shape_string = graph.is_completed(key)
    ? std::string("shape: ") + std::visit([] (auto p) { return to_string(p); }, graph.product_at(key))
    : std::string("");

    auto error_string = graph.is_error(key)
    ? std::string("error: ") + graph.error_at(key)
    : std::string("");

    os
    << std::left
    << std::setw(36)
    << std::setfill('.')
    << to_string(key)
    << " status: "
    << to_string(graph.status(key))
    << ' '
    << shape_string
    << error_string;
}

std::string represent_graph_item(const DependencyGraph& graph, product_identifier_t key)
{
    std::stringstream ss;
    represent_graph_item(ss, graph, key);
    return ss.str();
}

std::vector<std::string> represent_graph_items(const DependencyGraph& graph)
{
    return seq::view(graph.keys())
    | seq::map([&graph] (auto key) { return represent_graph_item(graph, key); })
    | seq::to<std::vector>();
}




//=============================================================================
auto basic_primitive(position_t p, mhd::magnetic_field_vector_t b)
{
    return mhd::primitive(1.0, {}, 1.0, b);
}

auto abc_vector_potential(position_t p)
{
    auto k = 2.0 * M_PI / dimensional::unit_length(1.0);
    auto [A, B, C] = std::tuple(1.0, 1.0, 1.0);
    auto [x, y, z] = as_tuple(p);
    auto b0 = 1.0;
    auto ax = A * std::sin(k * z) + C * std::cos(k * y);
    auto ay = B * std::sin(k * x) + A * std::cos(k * z);
    auto az = C * std::sin(k * y) + B * std::cos(k * x);
    return b0 * mhd::vector_potential_t{ax, ay, az};
};




//=============================================================================
auto build_graph()
{
    auto graph = DependencyGraph();
    auto depth = 2;
    auto block_size = 16;

    for (auto rule : mhd_rules::initial_condition_rules     (depth, block_size, basic_primitive, abc_vector_potential)) graph.define(rule);
    for (auto rule : mhd_rules::recover_primitive_rules     (depth, block_size, 0)) graph.define(rule);   
    for (auto rule : mhd_rules::primitive_extension_rules   (depth, block_size, 0)) graph.define(rule);
    for (auto rule : mhd_rules::magnetic_extension_rules    (depth, block_size, 0)) graph.define(rule);
    for (auto rule : mhd_rules::electromotive_force_rules   (depth, block_size, 0)) graph.define(rule);
    for (auto rule : mhd_rules::godunov_data_rules          (depth, block_size, 0)) graph.define(rule);
    for (auto rule : mhd_rules::global_primitive_array_rules(depth, block_size, 0)) graph.define(rule);

    return std::move(graph).throw_if_lacking_definitions();
}




//=============================================================================
int main()
{
    auto thread_count = 1;
    auto scheduler = mara::ThreadPool(thread_count);
    auto graph = build_graph();

    auto tb = ui::session_t();
    auto ui_state = ui::state_t();


    ui_state.content_table_size = [&graph] () { return graph.size(); };
    ui_state.content_table_item = [&graph, keys=graph.keys()] (unsigned row) { return represent_graph_item(graph, keys.at(row)); };
    ui_state.concurrent_task_count = [&scheduler] () { return scheduler.job_count(); };
    ui::draw(ui_state);


    while (true)
    {
        if (ui::fulfill(ui::action::quit))
        {
            break;
        }

        if (ui::fulfill(ui::action::evaluation_step))
        {
            for (const auto& key : graph.eligible_rules())
            {
                graph.evaluate_rule(key, scheduler);
            }
            ui::draw(ui_state);
        }

        if (ui::fulfill(ui::action::reset_simulation))
        {
            scheduler.restart(thread_count);
            graph = build_graph();
            ui::draw(ui_state);
        }

        if (! graph.poll(std::chrono::milliseconds(0)).empty())
        {
            ui::draw(ui_state);
        }

        if (auto event = ui::poll(std::chrono::milliseconds(5)); event.has_value())
        {
            ui::draw(ui_state = ui::handle_event(ui_state, event.value()));
        }
    }
    return 0;
}
