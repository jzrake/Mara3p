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
// #include "app_hdf5.hpp"
// #include "app_hdf5_dimensional.hpp"
// #include "app_hdf5_geometric.hpp"
// #include "app_hdf5_ndarray.hpp"
// #include "app_hdf5_numeric_array.hpp"
// #include "app_hdf5_ndarray_dimensional.hpp"
// #include "app_hdf5_rational.hpp"
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
void represent_graph_item(std::ostream& os, const DependencyGraph& graph, unsigned row)
{
    auto key = graph.keys().at(row);

    auto shape_string = graph.is_completed(key)
    ? std::visit([] (auto p) { return to_string(p); }, graph.product_at(key))
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
    << ' ' << shape_string
    << ' ' << error_string;
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




static auto depth = 2;
static auto block_size = 32;




//=============================================================================
auto block_indexes(nd::uint depth)
{
    return seq::adapt(nd::index_space(mesh::block_extent(depth)))
    | seq::map([depth] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; });
}

void enter_initial_conditions_rules(DependencyGraph& graph)
{
    using namespace mhd_rules;

    for (auto index : block_indexes(depth))
        for (auto rule : initial_condition(index, block_size, basic_primitive, abc_vector_potential))
            graph.define(rule);

}

void enter_iteration_rules(DependencyGraph& graph, rational::number_t iteration)
{
    using namespace mhd_rules;

    auto dt = dimensional::unit_time(0.01);

    for (auto index : block_indexes(depth))
    {
        graph.define(recover_primitive   (index, iteration));
        graph.define(extend_primitive    (index, iteration, block_size));
        graph.define(extend_magnetic     (index, iteration, block_size));
        graph.define(electromotive_force (index, iteration));
        graph.define(godunov_data        (index, iteration));
        graph.define(update_conserved    (index, iteration + 1, dt));
        graph.define(update_magnetic     (index, iteration + 1, dt));
    }
}




//=============================================================================
class RuntimeCoordinator
{
public:
    void reset(DependencyGraph& graph)
    {
        graph.clear();
        iteration = 0;
    }

    void update_definitions(DependencyGraph& graph)
    {
        if (graph.empty())
        {
            enter_initial_conditions_rules(graph);
        }
        else
        {
            enter_iteration_rules(graph, iteration);
            iteration = iteration + 1;
        }
    }

    rational::number_t current_iteration() const
    {
        return iteration;
    }

private:
    rational::number_t iteration;
};




//=============================================================================
auto vtk_output_side_effects()
{
    auto keys = block_indexes(depth) | seq::map([] (auto index)
    {
        return key(data_field::cell_primitive_variables).block(index).iteration(10).id;
    });

    return block_indexes(depth)
    | seq::map([] (auto index)
    {
        return [index] (const product_t& product)
        {
            auto order = dot(mesh::to_uivec(index.coordinates), nd::strides_row_major(mesh::block_extent(depth)));
            auto fname = "primitive." + std::to_string(order) + ".vtk";
            auto out   = std::fstream(fname, std::fstream::out | std::fstream::trunc | std::fstream::binary);
            auto xv    = mesh::construct_vertices<dimensional::unit_length>(index.level, index.coordinates, block_size);
            auto bc    = std::get<cell_primitive_variables_t>(product) | nd::map(mhd::magnetic_field_vector) | nd::to_shared();

            vtk::write(out, "primitive_fields", nd::to_shared(xv), std::pair("B", bc));
        };
    })
    | seq::keys(keys);
}




//=============================================================================
class SideEffects
{
public:
    void apply(const std::map<product_identifier_t, product_t>& updated_items) const
    {
        auto side_effects = vtk_output_side_effects() | seq::to_dict<std::map>();

        for (const auto& [k, p] : updated_items)
        {
            if (side_effects.count(k))
            {
                side_effects.at(k)(p);
            }
        }
    }
};




//=============================================================================
void drive(mara::ThreadPool& thread_pool, DependencyGraph& graph, RuntimeCoordinator& coordinator, SideEffects& side_effects)
{
    if (! graph.count_unevaluated())
    {
        coordinator.update_definitions(graph);
        graph.collect_garbage();
    }
    else
    {
        for (const auto& key : graph.eligible_rules())
        {
            graph.evaluate_rule(key, thread_pool);
        }
    }
    side_effects.apply(graph.poll(std::chrono::milliseconds(0)));
}




//=============================================================================
int main()
{
    auto thread_count = 1;
    auto thread_pool  = mara::ThreadPool(thread_count);
    auto coordinator  = RuntimeCoordinator();
    auto graph        = DependencyGraph();
    auto side_effects = SideEffects();


    auto termbox      = ui::session_t();
    auto ui_state     = ui::state_t();

    ui_state.content_table_size    = [&graph] () { return graph.size(); };
    ui_state.content_table_item    = [&graph] (unsigned row) { return represent_graph_item(graph, row); };
    ui_state.concurrent_task_count = [&thread_pool] () { return thread_pool.job_count(); };
    ui::draw(ui_state);


    while (long(coordinator.current_iteration()) < 100)
    {       
        if (ui::fulfill(ui::action::evaluation_step))
        {
            drive(thread_pool, graph, coordinator, side_effects);
        }

        if (ui::fulfill(ui::action::reset_simulation))
        {
            thread_pool.reset(thread_count);
            coordinator.reset(graph);
            ui::draw(ui_state);
        }

        if (ui::fulfill(ui::action::quit))
        {
            break;
        }

        if (auto event = ui::poll(std::chrono::milliseconds(20)); event.has_value())
        {
            ui_state = ui::handle_event(ui_state, event.value());
        }
        ui::draw(ui_state);
    }

    return 0;
}
