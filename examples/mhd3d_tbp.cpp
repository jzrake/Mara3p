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




#include <fstream>
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
#include "app_vtk.hpp"
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
void represent_graph_item(std::ostream& os, const DependencyGraph& graph, unsigned row)
{
    auto key = graph.keys().at(row);

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

std::string represent_graph_item(const DependencyGraph& graph, unsigned row)
{
    std::stringstream ss;
    represent_graph_item(ss, graph, row);
    return ss.str();
}

std::vector<std::string> represent_graph_items(const DependencyGraph& graph)
{
    return seq::range(graph.size())
    | seq::map([&graph] (auto row) { return represent_graph_item(graph, row); })
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




static auto depth = 1;
static auto block_size = 32;




//=============================================================================
void enter_initial_conditions_rules(DependencyGraph& graph)
{
    using namespace mhd_rules;

    for (auto rule : initial_condition_rules(depth, block_size, basic_primitive, abc_vector_potential))
    {
        graph.define(rule);
    }
}

void enter_iteration_rules(DependencyGraph& graph, rational::number_t iteration)
{
    using namespace mhd_rules;

    auto dt = dimensional::unit_time(0.1);

    for (auto rule : recover_primitive_rules     (depth, block_size, iteration)) graph.define(rule);   
    for (auto rule : primitive_extension_rules   (depth, block_size, iteration)) graph.define(rule);
    for (auto rule : magnetic_extension_rules    (depth, block_size, iteration)) graph.define(rule);
    for (auto rule : electromotive_force_rules   (depth, block_size, iteration)) graph.define(rule);
    for (auto rule : godunov_data_rules          (depth, block_size, iteration)) graph.define(rule);
    for (auto rule : update_conserved_rules      (depth, block_size, iteration + 1, dt)) graph.define(rule);
    for (auto rule : update_magnetic_rules       (depth, block_size, iteration + 1, dt)) graph.define(rule);
}

void reset(DependencyGraph& graph)
{
    graph = DependencyGraph();
    enter_initial_conditions_rules(graph);
}




//=============================================================================
void step_graph(DependencyGraph& graph, mara::ThreadPool& thread_pool)
{
    static int iteration = 0;

    if (graph.eligible_rules().empty() && iteration < 3)
    {
        enter_iteration_rules(graph, iteration);
        iteration += 1;
    }
    else
    {
        for (const auto& key : graph.eligible_rules())
        {
            graph.evaluate_rule(key, thread_pool);
        }
    }
}




//=============================================================================
auto vtk_output_side_effects()
{
    auto keys = mhd_rules::block_indexes(depth, block_size)
    | seq::map([] (auto index)
    {
        return key(data_field::cell_primitive_variables).block(index).id;
    });

    return block_indexes(depth, block_size)
    | seq::map([] (auto index)
    {
        return [index] (const product_t& product)
        {
            auto order = dot(mesh::to_uivec(index.coordinates), nd::strides_row_major(mesh::block_extent(depth)));
            auto fname = "prim." + std::to_string(order) + ".vtk";
            auto out   = std::fstream(fname, std::fstream::out | std::fstream::trunc | std::fstream::binary);
            auto xv    = mesh::construct_vertices<dimensional::unit_length>(index.level, index.coordinates, block_size);
            auto bc    = std::get<cell_primitive_variables_t>(product) | nd::map(mhd::magnetic_field_vector) | nd::to_shared();

            vtk::write(out, "primitive_fields", nd::to_shared(xv), std::pair("B", bc));
        };
    })
    | seq::keys(keys);
}




//=============================================================================
template<typename FunctionType>
bool side_effects(const std::map<product_identifier_t, FunctionType>& side_effects, const std::map<product_identifier_t, product_t>& updated_items)
{
    for (const auto& [k, p] : updated_items)
    {
        if (side_effects.count(k))
        {
            side_effects.at(k)(p);
        }
    }
    return ! updated_items.empty();
}




//=============================================================================
int main()
{
    auto thread_count = 12;
    auto thread_pool = mara::ThreadPool(thread_count);
    auto graph = DependencyGraph();


    auto side_effects_map = vtk_output_side_effects()
    | seq::to_dict<std::map>();


    auto tb = ui::session_t();
    auto ui_state = ui::state_t();

    reset(graph);

    ui_state.content_table_size = [&graph] () { return graph.size(); };
    ui_state.content_table_item = [&graph] (unsigned row) { return represent_graph_item(graph, row); };
    ui_state.concurrent_task_count = [&thread_pool] () { return thread_pool.job_count(); };
    ui::draw(ui_state);

    while (true)
    {
        if (ui::fulfill(ui::action::quit))
        {
            break;
        }

        if (ui::fulfill(ui::action::reset_simulation))
        {
            thread_pool.restart(thread_count);
            reset(graph);
            ui::draw(ui_state);
        }

        if (ui::fulfill(ui::action::evaluation_step))
        {
            step_graph(graph, thread_pool);
            ui::draw(ui_state);
        }

        if (side_effects(side_effects_map, graph.poll(std::chrono::milliseconds(0))))
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
