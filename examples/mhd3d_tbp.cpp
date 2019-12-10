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




#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <variant>
// #include "app_hdf5.hpp"
// #include "app_hdf5_dimensional.hpp"
// #include "app_hdf5_geometric.hpp"
// #include "app_hdf5_ndarray.hpp"
// #include "app_hdf5_numeric_array.hpp"
// #include "app_hdf5_ndarray_dimensional.hpp"
// #include "app_hdf5_rational.hpp"
#include "app_ui.hpp"
#include "core_bqo_tree.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_rational.hpp"
#include "core_sequence.hpp"
#include "mesh_cartesian_3d.hpp"
#include "parallel_dependency_graph.hpp"
#include "parallel_thread_pool.hpp"
#include "physics_mhd.hpp"
#include "scheme_mhd_rules.hpp"
#include "scheme_mhd_v2.hpp"




using namespace mhd_scheme_v2;
using namespace mhd_rules;




//=============================================================================
namespace h5 {

// template<typename T>
// void write(const Group& group, std::string name, const std::array<nd::shared_array<T, 3>, 3>& v)
// {
//     write(group.require_group(name), "1", v.at(0));
//     write(group.require_group(name), "2", v.at(1));
//     write(group.require_group(name), "3", v.at(2));
// }

// void write(const Group& group, std::string name, const product_t& product)
// {
//     std::visit([&group, name] (auto p) { write(group, name, p); }, product);
// }

// std::string legalize(std::string s)
// {
//     return seq::view(s)
//     | seq::map([] (auto c) { return c == '/' ? '@' : c; })
//     | seq::remove_if([] (auto c) { return c == ' '; })
//     | seq::to<std::basic_string>();
// }

}




//=============================================================================
std::string to_string(mara::evaluation_status s)
{
    switch (s)
    {
        case mara::evaluation_status::undefined:   return "?";
        case mara::evaluation_status::error:       return "!";
        case mara::evaluation_status::defined:     return "d";
        case mara::evaluation_status::pending:     return ".";
        case mara::evaluation_status::eligible:    return "e";
        case mara::evaluation_status::completed:   return "c";
        default: return "";
    }
}

std::string to_string(rational::number_t n)
{
    return std::to_string(n.num) + "/" + std::to_string(n.den);
}

std::string to_string(multilevel_index_t i)
{
    return std::to_string(i.level)
    + ":["
    + std::to_string(i.coordinates[0]) + ", "
    + std::to_string(i.coordinates[1]) + ", "
    + std::to_string(i.coordinates[2]) + "]";
}

std::string to_string(data_field f)
{
    switch (f)
    {
        case data_field::cell_primitive_variables:     return "P";
        case data_field::cell_conserved_density:       return "U";
        case data_field::face_godunov_data:            return "F";
        case data_field::face_magnetic_flux_density:   return "B";
        case data_field::edge_electromotive_density:   return "E";
        case data_field::edge_vector_potential:        return "A";
        default: return "";
    }
}

std::string to_string(extended_status s)
{
    switch (s)
    {
        case extended_status::not_extended:  return " [] ";
        case extended_status::extended:      return "[[]]";
        default: return "";
    }
}

std::string to_string(product_identifier_t id)
{
    auto [a, b, c, d] = id;
    return
      to_string(a) + " - "
    + to_string(b) + " - "
    + to_string(c) + " - "
    + to_string(d);
}




//=============================================================================
template<typename T>
auto represent_shape(const nd::shared_array<T, 3>& v)
{
    return "("
    + std::to_string(shape(v, 0)) + " "
    + std::to_string(shape(v, 1)) + " "
    + std::to_string(shape(v, 2)) + ")";
}

template<typename T>
auto represent_shape(const std::array<nd::shared_array<T, 3>, 3>& v)
{
    return represent_shape(v.at(0)) + " " + represent_shape(v.at(1)) + " " + represent_shape(v.at(2));
}

void represent_graph_item(std::ostream& os, const DependencyGraph& graph, product_identifier_t key)
{
    auto shape_string = graph.is_completed(key)
    ? std::string("shape: ") + std::visit([] (auto p) { return represent_shape(p); }, graph.product_at(key))
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

std::vector<std::string> represent_graph_items(const DependencyGraph& graph)
{
    return seq::view(graph.keys())
    | seq::map([&graph] (auto key)
    {
        std::stringstream ss;
        represent_graph_item(ss, graph, key);
        return ss.str();
    })
    | seq::to<std::vector>();
}




//=============================================================================
auto build_graph()
{
    auto graph = DependencyGraph();
    auto depth = 2;
    auto block_size = 16;

    for (auto rule : mhd_rules::initial_condition_rules     (depth, block_size))    graph.define(rule);
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

    ui_state.content_table_items = represent_graph_items(graph);
    ui_state.concurrent_task_count = scheduler.job_count();
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
        }

        if (ui::fulfill(ui::action::reset_simulation))
        {
            scheduler.restart(thread_count);
            graph = build_graph();

            ui_state.content_table_items = represent_graph_items(graph);
            ui_state.concurrent_task_count = scheduler.job_count();
            ui::draw(ui_state);
        }

        if (! graph.poll(std::chrono::milliseconds(5)).empty())
        {
            ui_state.content_table_items = represent_graph_items(graph);
            ui_state.concurrent_task_count = scheduler.job_count();
            ui::draw(ui_state);
        }

        if (auto event = ui::peek(5); event.has_value())
        {
            ui::draw(ui_state = ui::handle_event(ui_state, event.value()));
        }
    }

    return 0;
}




//=============================================================================
// void execute(DependencyGraph& graph)
// {
//     auto scheduler = mara::ThreadPool(1);
//     auto start = std::chrono::high_resolution_clock::now();

//     while (graph.count_unevaluated(is_responsible_for))
//     {
//         for (const auto& key : graph.eligible_rules(is_responsible_for))
//         {
//             graph.evaluate_rule(key, scheduler);
//         }

//         if (! graph.poll(std::chrono::milliseconds(50)).empty())
//         {
//             print_graph_status(graph, is_responsible_for);
//         }
//     }

//     auto delta = std::chrono::high_resolution_clock::now() - start;

//     std::cout
//     << "Time to complete graph: "
//     << 1e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(delta).count()
//     << "s"
//     << std::endl;

//     auto h5f = h5::File("test.h5", "w");

//     for (const auto& [key, product] : graph.items())
//     {
//         h5::write(h5f, h5::legalize(to_string(key)), product);

//         if (std::get<1>(key).level == 0)
//         {
//             h5::write(h5f, "primitive", product);
//         }
//     }
// }
