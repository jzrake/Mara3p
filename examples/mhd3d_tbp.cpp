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
#include "scheme_mhd_v2.hpp"




using namespace mhd_scheme_v2;




//=============================================================================
enum class data_field
{
    cell_primitive_variables,
    cell_conserved_density,
    face_godunov_data,
    face_magnetic_flux_density,
    edge_electromotive_density,
    edge_vector_potential,
};

enum class extended_status
{
    not_extended,
    extended,
};

using product_t = std::variant<
    cell_primitive_variables_t,
    cell_conserved_density_t,
    face_godunov_data_t,
    face_magnetic_flux_density_t,
    edge_electromotive_density_t>;

using product_identifier_t = std::tuple<
    rational::number_t,
    multilevel_index_t,
    data_field,
    extended_status>;

using DependencyGraph = mara::DependencyGraph<product_identifier_t, product_t>;




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
        case mara::evaluation_status::undefined:   return "!";
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
    auto shape_string = graph.is_completed(key) ? std::visit([] (auto p) { return represent_shape(p); }, graph.product_at(key)) : std::string("?");

    os
    << std::left
    << std::setw(36)
    << std::setfill('.')
    << to_string(key)
    << " status: "
    << to_string(graph.status(key))
    << " shape: "
    << shape_string;
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
static auto is_responsible_for = [] (auto) { return true; };
static auto block_size = 16;
static auto depth = 1UL;
static auto block_extent = nd::uivec(1 << depth, 1 << depth, 1 << depth);




//=============================================================================
using named_rule_t = std::tuple<
    product_identifier_t,
    DependencyGraph::mapping_type,
    std::vector<product_identifier_t>>;

named_rule_t rule(
    product_identifier_t key,
    DependencyGraph::mapping_type mapping,
    std::vector<product_identifier_t> args)
{
    return named_rule_t{key, mapping, args};
}




//=============================================================================
struct key_factory_t
{
    key_factory_t iteration(rational::number_t v)                        const { auto n = id; std::get<0>(n) = v; return {n}; }
    key_factory_t block    (multilevel_index_t v)                        const { auto n = id; std::get<1>(n) = v; return {n}; }
    key_factory_t field    (data_field v)                                const { auto n = id; std::get<2>(n) = v; return {n}; }
    key_factory_t extended (extended_status v=extended_status::extended) const { auto n = id; std::get<3>(n) = v; return {n}; }

    operator product_identifier_t() const
    {
        return id;
    }

    auto bind_block() const { return [this] (multilevel_index_t v) { return block(v).id; }; }

    product_identifier_t id;
};

key_factory_t key(data_field field)
{
    return key_factory_t().field(field);
}

key_factory_t key(rational::number_t iteration)
{
    return key_factory_t().iteration(iteration);
}




//=============================================================================
template<typename FunctionType>
auto wrap(FunctionType f)
{
    return [f] (std::vector<product_t>) -> product_t
    {
        return f();
    };
}

template<typename Arg1, typename FunctionType>
auto wrap(FunctionType f)
{
    return [f] (std::vector<product_t> args) -> product_t
    {
        const auto& arg1 = std::get<Arg1>(args.at(0));
        return f(arg1);
    };
}

template<typename Arg1, typename Arg2, typename FunctionType>
auto wrap(FunctionType f)
{
    return [f] (std::vector<product_t> args) -> product_t
    {
        const auto& arg1 = std::get<Arg1>(args.at(0));
        const auto& arg2 = std::get<Arg2>(args.at(1));
        return f(arg1, arg2);
    };
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
auto block_indexes()
{
    return seq::adapt(nd::index_space(block_extent))
    | seq::map([] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; });
}

auto extend_cell_primitive_variables(nd::uint block_size, nd::uint count)
{
    return [=] (const std::vector<product_t>& vargs) -> product_t
    {
        auto pargs = seq::get<cell_primitive_variables_t>(seq::view(vargs)) | seq::to<std::vector>();
        return extend_cell_primitive_variables(pargs, block_size, count);
    };
}

auto extend_face_magnetic_flux_density(nd::uint block_size, nd::uint count)
{
    return [=] (const std::vector<product_t>& vargs) -> product_t
    {
        auto pargs = seq::get<face_magnetic_flux_density_t>(seq::view(vargs)) | seq::to<std::vector>();
        return extend_face_magnetic_flux_density(pargs, block_size, count);
    };
}




//=============================================================================
auto initial_condition_rules()
{
    auto nb = unsigned(1 << depth);
    auto dl = dimensional::unit_length(1.0 / block_size / nb);

    return block_indexes()
    | seq::map([dl] (auto index)
    {
        auto k_ae = key(data_field::edge_vector_potential)     .block(index);
        auto k_bf = key(data_field::face_magnetic_flux_density).block(index);
        auto k_uc = key(data_field::cell_conserved_density)    .block(index);

        auto f_ae = wrap([index] () { return construct_vector_potential(index, block_size, abc_vector_potential); });
        auto f_bf = wrap<edge_electromotive_density_t>([dl] (auto ee) { return curl(ee, dl); });
        auto f_uc = wrap<face_magnetic_flux_density_t>([index] (auto bf) { return construct_conserved(index, block_size, bf, basic_primitive); });

        auto r_ae = rule(k_ae, f_ae, {});
        auto r_bf = rule(k_bf, f_bf, {k_ae});
        auto r_uc = rule(k_uc, f_uc, {k_bf});

        return seq::from(r_ae, r_bf, r_uc);
    })
    | seq::flat()
    | seq::to_dynamic();
}




//=============================================================================
auto recover_primitive_rules(rational::number_t iteration)
{
    return block_indexes()
    | seq::map([iteration] (auto index)
    {
        auto bf = key(iteration).block(index).field(data_field::face_magnetic_flux_density);
        auto uc = key(iteration).block(index).field(data_field::cell_conserved_density);
        auto pc = key(iteration).block(index).field(data_field::cell_primitive_variables);
        auto fm = wrap<cell_conserved_density_t, face_magnetic_flux_density_t>(mhd_scheme_v2::primitive_array);
        return rule(pc, fm, {uc, bf});
    })
    | seq::to_dynamic();
}




//=============================================================================
auto primitive_extension_rules(rational::number_t iteration)
{
    return block_indexes()
    | seq::map([iteration] (auto index) -> named_rule_t
    {
        auto pc = key(data_field::cell_primitive_variables).iteration(iteration);

        auto neighbor_ids = mesh::neighbors_27(mesh::to_uivec(index.coordinates), block_extent)
        | seq::map([] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; })
        | seq::map(pc.bind_block())
        | seq::to<std::vector>();

        return rule(pc.extended().block(index), extend_cell_primitive_variables(block_size, 1), neighbor_ids);
    })
    | seq::to_dynamic();
}




//=============================================================================
auto magnetic_extension_rules(rational::number_t iteration)
{
    return block_indexes()
    | seq::map([iteration] (auto index) -> named_rule_t
    {
        auto bf = key(data_field::face_magnetic_flux_density).iteration(iteration);
        auto n1 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), block_extent, mesh::axis_3d::i);
        auto n2 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), block_extent, mesh::axis_3d::j);
        auto n3 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), block_extent, mesh::axis_3d::k);

        auto neighbor_ids = seq::concat(n1, n2, n3)
        | seq::map([] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; })
        | seq::map(bf.bind_block())
        | seq::to<std::vector>();

        return rule(bf.extended().block(index), extend_face_magnetic_flux_density(block_size, 1), neighbor_ids);
    })
    | seq::to_dynamic();
}




//=============================================================================
auto godunov_data_rules(rational::number_t iteration)
{
    return block_indexes()
    | seq::map([iteration] (auto index) -> named_rule_t
    {
        auto pc = key(data_field::cell_primitive_variables)  .iteration(iteration).block(index).extended();
        auto bf = key(data_field::face_magnetic_flux_density).iteration(iteration).block(index).extended();
        auto gf = key(data_field::face_godunov_data)         .iteration(iteration).block(index);
        auto fm = wrap<cell_primitive_variables_t, face_magnetic_flux_density_t>(godunov_fluxes);
        return rule(gf, fm, {pc, bf});
    })
    | seq::to_dynamic();
}




//=============================================================================
auto electromotive_force_rules(rational::number_t iteration)
{
    return block_indexes()
    | seq::map([iteration] (auto index) -> named_rule_t
    {
        auto ee = key(data_field::edge_electromotive_density).iteration(iteration).block(index);
        auto gf = key(data_field::face_godunov_data         ).iteration(iteration).block(index);
        auto fm = wrap<face_godunov_data_t>(electromotive_forces);
        return rule(ee, fm, {gf});
    })
    | seq::to_dynamic();
}




//=============================================================================
auto global_primitive_array_rules(rational::number_t iteration)
{
    auto pc = key(data_field::cell_primitive_variables).iteration(iteration);
    auto arg_keys = block_indexes() | seq::map(pc.bind_block()) | seq::to<std::vector>();

    auto tile = [] (std::vector<product_t> block_vector) -> product_t
    {
        auto block_map = seq::view(block_vector)
        | seq::map([] (auto p) { return std::get<cell_primitive_variables_t>(p); })
        | seq::keys(nd::index_space(block_extent))
        | seq::to_dict<std::map>();

        return mesh::tile_blocks(block_map, block_extent) | nd::to_shared();
    };
    return seq::just(rule(pc, tile, arg_keys)) | seq::to_dynamic();
}




//=============================================================================
auto build_graph()
{
    auto graph = DependencyGraph();

    for (auto rule : initial_condition_rules())       graph.define(rule);
    for (auto rule : recover_primitive_rules(0))      graph.define(rule);   
    for (auto rule : primitive_extension_rules(0))    graph.define(rule);
    for (auto rule : magnetic_extension_rules(0))     graph.define(rule);
    for (auto rule : electromotive_force_rules(0))    graph.define(rule);
    for (auto rule : godunov_data_rules(0))           graph.define(rule);
    for (auto rule : global_primitive_array_rules(0)) graph.define(rule);

    return std::move(graph).throw_if_incomplete();
}




//=============================================================================
int main()
{
    auto scheduler = mara::ThreadPool(1);
    auto graph = build_graph();

    auto tb = ui::session_t();
    auto ui_state = ui::state_t();

    ui_state.content_table_items = represent_graph_items(graph);
    ui::draw(ui_state);


    while (true)
    {
        if (ui::fulfill(ui::action::quit))
        {
            break;
        }

        if (ui::fulfill(ui::action::evaluation_step))
        {
            for (const auto& key : graph.eligible_rules(is_responsible_for))
            {
                graph.evaluate_rule(key, scheduler);
            }
        }

        if (ui::fulfill(ui::action::reset_simulation))
        {
            scheduler.restart(1);
            graph = build_graph();

            ui_state.content_table_items = represent_graph_items(graph);
            ui::draw(ui_state);
        }

        if (! graph.poll(std::chrono::milliseconds(5)).empty())
        {
            ui_state.content_table_items = represent_graph_items(graph);
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

