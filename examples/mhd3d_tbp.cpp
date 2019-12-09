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
#include "app_hdf5.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_geometric.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_rational.hpp"
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
#include "scheme_mhd.hpp"



#if 1

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




//=============================================================================
using cell_primitive_variables_t   = nd::shared_array<mhd::primitive_t, 3>;
using cell_conserved_density_t     = nd::shared_array<mhd::conserved_density_t, 3>;
using face_godunov_data_t          = std::array<nd::shared_array<mhd::godunov_data_t, 3>, 3>;
using face_magnetic_flux_density_t = std::array<nd::shared_array<mhd::unit_magnetic_field, 3>, 3>;
using edge_electromotive_density_t = std::array<nd::shared_array<mhd::unit_electric_field, 3>, 3>;
using position_t                   = geometric::euclidean_vector_t<dimensional::unit_length>;
using multilevel_index_t           = bqo_tree::tree_index_t<3>;

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

template<typename T>
void write(const Group& group, std::string name, const std::array<nd::shared_array<T, 3>, 3>& v)
{
    write(group.require_group(name), "1", v.at(0));
    write(group.require_group(name), "2", v.at(1));
    write(group.require_group(name), "3", v.at(2));
}

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
static auto is_responsible_for = [] (auto) { return true; };
static auto block_size = 16;
static auto depth = 1UL;
static auto block_extent = nd::uivec(1 << depth, 1 << depth, 1 << depth);
static auto gamma_law_index = 5. / 3;




using namespace std::placeholders;
using geometric::unit_vector_on;
using util::apply_to;




// These are riemann solver functions for each of the three axes.
static auto rs1 = apply_to(std::bind(mhd::riemann_hlle, _1, _2, _3, unit_vector_on(1), gamma_law_index));
static auto rs2 = apply_to(std::bind(mhd::riemann_hlle, _1, _2, _3, unit_vector_on(2), gamma_law_index));
static auto rs3 = apply_to(std::bind(mhd::riemann_hlle, _1, _2, _3, unit_vector_on(3), gamma_law_index));




// These are "remove-transverse" functions; they return a rectangular pipe from
// the middle of a 3D array along a given axis.
static auto rt1 = mesh::remove_transverse_i(1);
static auto rt2 = mesh::remove_transverse_j(1);
static auto rt3 = mesh::remove_transverse_k(1);




// These are "remove-longitudinal" functions; they remove a layer of zones from
// the ends of a 3D array along a given axis.
static auto rl1 = nd::select(0, 1, -1);
static auto rl2 = nd::select(1, 1, -1);
static auto rl3 = nd::select(2, 1, -1);




//=============================================================================
void graph_item(std::ostream& os, const DependencyGraph& graph, product_identifier_t key)
{
    os
    << std::left
    << std::setw(48)
    << std::setfill('.')
    << to_string(key)
    << " status: "
    << to_string(graph.status(key))
    << '\n';
}

std::vector<std::string> represent_graph_items(const DependencyGraph& graph)
{
    return seq::view(graph.keys())
    | seq::map([&graph] (auto key)
    {
        std::stringstream ss;
        graph_item(ss, graph, key);
        return ss.str();
    })
    | seq::to<std::vector>();
}




//=========================================================================
template<typename IsResponsibleFor>
void print_graph_status(const DependencyGraph& graph, IsResponsibleFor is_responsible_for)
{
    std::invoke([&] ()
    {
        std::cout
        << std::string(52, '=')
        << "\nGraph status ("
        << graph.count_unevaluated(is_responsible_for)
        << " unevaluated):\n\n";

        for (auto key : graph.keys())
        {
            graph_item(std::cout, graph, key);
        }
        std::cout << '\n';
    });
}




//=============================================================================
using named_rule_t = std::tuple<
    product_identifier_t,
    std::function<product_t(std::vector<product_t>)>,
    std::vector<product_identifier_t>>;

named_rule_t rule(
    product_identifier_t key,
    std::function<product_t(std::vector<product_t>)> mapping,
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
auto construct_vertices(multilevel_index_t index)
{
    auto N = block_size;
    auto dx = dimensional::unit_length(1.0 / (1 << index.level));
    auto x0 = dx * geometric::to_euclidean_vector(numeric::construct<double>(index.coordinates));
    auto xv = dx * mesh::unit_lattice(N + 1, N + 1, N + 1) + x0;
    return xv;
}




//=============================================================================
auto construct_vector_potential(multilevel_index_t index)
{
    return [index] ()
    {
        auto A = [] (unsigned dir) { return util::compose(geometric::component(dir), abc_vector_potential); };
        auto dt = dimensional::unit_time(1.0);
        auto [xe1, xe2, xe3] = mesh::edge_positions(construct_vertices(index));

        return edge_electromotive_density_t{
            ((xe1 | nd::map(A(1))) / dt) | nd::to_shared(),
            ((xe2 | nd::map(A(2))) / dt) | nd::to_shared(),
            ((xe3 | nd::map(A(3))) / dt) | nd::to_shared(),
        };
    };
}




//=============================================================================
auto construct_conserved(multilevel_index_t index)
{
    return [index] (face_magnetic_flux_density_t bf)
    {
        auto p2c = std::bind(mhd::conserved_density, std::placeholders::_1, gamma_law_index);
        auto [bf1, bf2, bf3] = bf;
        auto xc = mesh::cell_positions(construct_vertices(index));
        auto bc = mesh::face_to_cell(bf1, bf2, bf3);
        auto pc = nd::zip(xc, bc) | nd::map(util::apply_to(basic_primitive));
        auto uc = pc | nd::map(p2c);
        return uc | nd::to_shared();
    };
}




//=============================================================================
auto curl(dimensional::unit_length dl)
{
    return [dl] (edge_electromotive_density_t ee)
    {
        auto dt = dimensional::unit_time(1.0);
        auto [ee1, ee2, ee3] = ee;
        auto [cf1, cf2, cf3] = mesh::solenoidal_difference(ee1, ee2, ee3);

        return face_magnetic_flux_density_t{
            cf1 * dt / dl | nd::to_shared(),
            cf2 * dt / dl | nd::to_shared(),
            cf3 * dt / dl | nd::to_shared(),
        };
    };
}




//=============================================================================
auto extend_cell_primitive_variables(nd::uint block_size, nd::uint count)
{
    return [r=block_size - count] (std::vector<product_t> args) -> product_t
    {
        auto items = seq::view(args)
        | seq::map([] (auto p) { return std::get<cell_primitive_variables_t>(p); })
        | seq::keys(nd::index_space(3, 3, 3))
        | seq::to_dict<std::map>();

        return mesh::tile_blocks(items, nd::uivec(3, 3, 3)) | mesh::remove_surface(r) | nd::to_shared();
    };
}




//=============================================================================
auto extend_face_magnetic_flux_density(nd::uint block_size, nd::uint count)
{
    return [r=block_size - count] (std::vector<product_t> args) -> product_t
    {
        auto bfs = seq::get<face_magnetic_flux_density_t>(seq::view(args)) | seq::chunk(9);
        auto axs = seq::from(mesh::axis_3d::i, mesh::axis_3d::j, mesh::axis_3d::k);
        auto shp = seq::from(nd::uivec(1, 3, 3), nd::uivec(3, 1, 3), nd::uivec(3, 3, 1));

        auto extend = [r] (unsigned a, auto bf, mesh::axis_3d axis, nd::uivec_t<3> shape)
        {
            auto items = seq::adapt(bf)
            | seq::map([a] (auto x) { return x.at(a); })
            | seq::keys(nd::index_space(shape))
            | seq::to_dict<std::map>();

            return mesh::tile_blocks(items, shape)
            | mesh::remove_transverse(r, axis)
            | nd::to_shared();
        };

        return seq::zip(seq::range(3), bfs, axs, shp)
        | seq::map(util::apply_to(extend))
        | seq::to_std_array<3>();
    };
}




/**
 * @brief      Generate Godunov fluxes and EMF's using piecewise constant
 *             (first-order) extrapolation.
 *
 * @param[in]  pc    Cell-centered primitive quantities, extended by 1 zones on
 *                   all three axes
 * @param[in]  bf    Face B-fields, extended by 1 zone on the transverse axes
 *
 * @return     A 3-array of face-centered Godunov fluxes and EMF's. Each of the
 *             returned arrays are extended by 1 zone on the transverse axes.
 */
face_godunov_data_t godunov_fluxes(cell_primitive_variables_t pc, face_magnetic_flux_density_t bf)
{
    auto to_godunov_data = apply_to([] (auto ff, auto ef) -> mhd::godunov_data_t
    {
        return numeric::tuple_cat(ff, numeric::tuple(ef.component_1(), ef.component_2(), ef.component_3()));
    });

    auto [bf1, bf2, bf3] = bf;

    auto pf1 = nd::zip(pc | nd::select(0, 0, -1), pc | nd::select(0, 1), bf1);
    auto pf2 = nd::zip(pc | nd::select(1, 0, -1), pc | nd::select(1, 1), bf2);
    auto pf3 = nd::zip(pc | nd::select(2, 0, -1), pc | nd::select(2, 1), bf3);

    auto gf1 = pf1 | nd::map(rs1) | nd::map(to_godunov_data) | nd::to_shared();
    auto gf2 = pf2 | nd::map(rs2) | nd::map(to_godunov_data) | nd::to_shared();
    auto gf3 = pf3 | nd::map(rs3) | nd::map(to_godunov_data) | nd::to_shared();

    return std::array{gf1, gf2, gf3};
}




/**
 * @brief      Generate EMF values on cell edges by averaging the electric field
 *             estimates from godunov data on the cell faces.
 *
 * @param[in]  gf    The godunoc data on cell faces
 *
 * @return     A 3-array of electric field values at edge midpoints, in the
 *             direction of the respective edge.
 */
edge_electromotive_density_t electromotive_forces(face_godunov_data_t gf)
{
    auto c1 = nd::map([] (auto g) { return numeric::get<5>(g); });
    auto c2 = nd::map([] (auto g) { return numeric::get<6>(g); });
    auto c3 = nd::map([] (auto g) { return numeric::get<7>(g); });
    auto a1 = nd::adjacent_mean(0);
    auto a2 = nd::adjacent_mean(1);
    auto a3 = nd::adjacent_mean(2);
    auto ev = nd::to_shared();
    auto [gf1, gf2, gf3] = gf;

    auto ee1a = gf2 | c1 | a3;
    auto ee2a = gf3 | c2 | a1;
    auto ee3a = gf1 | c3 | a2;
    auto ee1b = gf3 | c1 | a2;
    auto ee2b = gf1 | c2 | a3;
    auto ee3b = gf2 | c3 | a1;

    auto ee1 = 0.5 * (ee1a + ee1b);
    auto ee2 = 0.5 * (ee2a + ee2b);
    auto ee3 = 0.5 * (ee3a + ee3b);

    return std::array{ee1|ev, ee2|ev, ee3|ev};
}




/**
 * @brief      { function_description }
 *
 * @param[in]  uc    { parameter_description }
 * @param[in]  gf    { parameter_description }
 *
 * @return     { description_of_the_return_value }
 */
cell_conserved_density_t updated_conserved_density(cell_conserved_density_t uc, face_godunov_data_t gf)
{
    auto conserved_density_components = [] (auto g)
    {
        return numeric::tuple(
            numeric::get<0>(g),
            numeric::get<1>(g),
            numeric::get<2>(g),
            numeric::get<3>(g),
            numeric::get<4>(g));
    };

    auto dt = dimensional::unit_time(1.0);
    auto dl = dimensional::unit_length(1.0);

    auto [gf1, gf2, gf3] = gf;
    auto ff1 = gf1 | nd::map(conserved_density_components) | rt1;
    auto ff2 = gf2 | nd::map(conserved_density_components) | rt2;
    auto ff3 = gf3 | nd::map(conserved_density_components) | rt3;
    return (uc - mesh::divergence_difference(ff1, ff2, ff3) * dt / dl) | nd::to_shared();
}




/**
 * @brief      { function_description }
 *
 * @param[in]  bf    { parameter_description }
 * @param[in]  ee    { parameter_description }
 *
 * @return     { description_of_the_return_value }
 */
face_magnetic_flux_density_t updated_magnetic_flux_density(face_magnetic_flux_density_t bf, edge_electromotive_density_t ee)
{
    auto [bf1, bf2, bf3] = bf;
    auto [ee1, ee2, ee3] = ee;
    auto [cf1, cf2, cf3] = mesh::solenoidal_difference(ee1|rl1, ee2|rl2, ee3|rl3);

    auto dt = dimensional::unit_time(1.0);
    auto dl = dimensional::unit_length(1.0);

    return std::array{
        (bf1 - cf1 * dt / dl) | nd::to_shared(),
        (bf2 - cf2 * dt / dl) | nd::to_shared(),
        (bf3 - cf3 * dt / dl) | nd::to_shared()};
}




//=============================================================================
auto block_indexes()
{
    return seq::adapt(nd::index_space(block_extent))
    | seq::map([] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; });
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

        auto f_ae = wrap(construct_vector_potential(index));
        auto f_bf = wrap<edge_electromotive_density_t>(curl(dl));
        auto f_uc = wrap<face_magnetic_flux_density_t>(construct_conserved(index));

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
        auto fm = wrap<cell_conserved_density_t, face_magnetic_flux_density_t>(mara::primitive_array);
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

    std::cout << "Building graph..." << std::endl;

    for (auto rule : initial_condition_rules())       graph.define(rule);
    for (auto rule : recover_primitive_rules(0))      graph.define(rule);   
    for (auto rule : primitive_extension_rules(0))    graph.define(rule);
    for (auto rule : magnetic_extension_rules(0))     graph.define(rule);
    for (auto rule : electromotive_force_rules(0))    graph.define(rule);
    for (auto rule : godunov_data_rules(0))           graph.define(rule);
    for (auto rule : global_primitive_array_rules(0)) graph.define(rule);

    std::cout << "Done" << std::endl;

    return std::move(graph).throw_if_incomplete();
}




//=============================================================================
void execute(DependencyGraph& graph)
{
    auto scheduler = mara::ThreadPool(1);
    auto start = std::chrono::high_resolution_clock::now();

    while (graph.count_unevaluated(is_responsible_for))
    {
        for (const auto& key : graph.eligible_rules(is_responsible_for))
        {
            graph.evaluate_rule(key, scheduler);
        }

        if (! graph.poll(std::chrono::milliseconds(50)).empty())
        {
            print_graph_status(graph, is_responsible_for);
        }
    }

    auto delta = std::chrono::high_resolution_clock::now() - start;

    std::cout
    << "Time to complete graph: "
    << 1e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(delta).count()
    << "s"
    << std::endl;

    auto h5f = h5::File("test.h5", "w");

    for (const auto& [key, product] : graph.items())
    {
        h5::write(h5f, h5::legalize(to_string(key)), product);

        if (std::get<1>(key).level == 0)
        {
            h5::write(h5f, "primitive", product);
        }
    }
}

#endif




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
        if (auto event = ui::peek(5); event.has_value())
        {
            if (ui::is_quit(event.value()))
            {
                break;
            }
            ui::draw(ui_state = ui::handle_event(ui_state, event.value()));
        }

        for (const auto& key : graph.eligible_rules(is_responsible_for))
        {
            graph.evaluate_rule(key, scheduler);
        }

        if (! graph.poll(std::chrono::milliseconds(5)).empty())
        {
            ui_state.content_table_items = represent_graph_items(graph);

            ui::draw(ui_state);
        }
    }

    return 0;
}
