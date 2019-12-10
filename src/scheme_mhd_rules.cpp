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




#include "scheme_mhd_rules.hpp"




using namespace mhd_scheme_v2;
using namespace mhd_rules;




//=============================================================================
static nd::uivec_t<3> block_extent(nd::uint depth)
{
    return nd::uivec(1 << depth, 1 << depth, 1 << depth);
}

static auto block_indexes(nd::uint depth, nd::uint block_size)
{
    return seq::adapt(nd::index_space(block_extent(depth)))
    | seq::map([depth] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; });
}

static auto extend_cell_primitive_variables(nd::uint block_size, nd::uint count)
{
    return [=] (const std::vector<product_t>& vargs) -> product_t
    {
        auto pargs = seq::get<cell_primitive_variables_t>(seq::view(vargs)) | seq::to<std::vector>();
        return extend_cell_primitive_variables(pargs, block_size, count);
    };
}

static auto extend_face_magnetic_flux_density(nd::uint block_size, nd::uint count)
{
    return [=] (const std::vector<product_t>& vargs) -> product_t
    {
        auto pargs = seq::get<face_magnetic_flux_density_t>(seq::view(vargs)) | seq::to<std::vector>();
        return extend_face_magnetic_flux_density(pargs, block_size, count);
    };
}




//=============================================================================
template<typename FunctionType>
static auto wrap(FunctionType f)
{
    return [f] (std::vector<product_t>) -> product_t
    {
        return f();
    };
}

template<typename Arg1, typename FunctionType>
static auto wrap(FunctionType f)
{
    return [f] (std::vector<product_t> args) -> product_t
    {
        const auto& arg1 = std::get<Arg1>(args.at(0));
        return f(arg1);
    };
}

template<typename Arg1, typename Arg2, typename FunctionType>
static auto wrap(FunctionType f)
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
seq::dynamic_sequence_t<named_rule_t> mhd_rules::initial_condition_rules(nd::uint depth, nd::uint block_size)
{
    auto nb = unsigned(1 << depth);
    auto dl = dimensional::unit_length(1.0 / block_size / nb);

    return block_indexes(depth, block_size)
    | seq::map([=] (auto index)
    {
        auto k_ae = key(data_field::edge_vector_potential)     .block(index);
        auto k_bf = key(data_field::face_magnetic_flux_density).block(index);
        auto k_uc = key(data_field::cell_conserved_density)    .block(index);

        auto f_ae = wrap([=] () { return construct_vector_potential(index, block_size, abc_vector_potential); });
        auto f_bf = wrap<edge_electromotive_density_t>([=] (auto ee) { return curl(ee, dl); });
        auto f_uc = wrap<face_magnetic_flux_density_t>([=] (auto bf) { return construct_conserved(index, block_size, bf, basic_primitive); });

        auto r_ae = rule(k_ae, f_ae, {});
        auto r_bf = rule(k_bf, f_bf, {k_ae});
        auto r_uc = rule(k_uc, f_uc, {k_bf});

        return seq::from(r_ae, r_bf, r_uc);
    })
    | seq::flat()
    | seq::to_dynamic();
}




//=============================================================================
seq::dynamic_sequence_t<named_rule_t> mhd_rules::recover_primitive_rules(nd::uint depth, nd::uint block_size, rational::number_t iteration)
{
    return block_indexes(depth, block_size)
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
seq::dynamic_sequence_t<named_rule_t> mhd_rules::primitive_extension_rules(nd::uint depth, nd::uint block_size, rational::number_t iteration)
{
    return block_indexes(depth, block_size)
    | seq::map([=] (auto index) -> named_rule_t
    {
        auto pc = key(data_field::cell_primitive_variables).iteration(iteration);

        auto neighbor_ids = mesh::neighbors_27(mesh::to_uivec(index.coordinates), block_extent(depth))
        | seq::map([=] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; })
        | seq::map(pc.bind_block())
        | seq::to<std::vector>();

        return rule(pc.extended().block(index), extend_cell_primitive_variables(block_size, 1), neighbor_ids);
    })
    | seq::to_dynamic();
}




//=============================================================================
seq::dynamic_sequence_t<named_rule_t> mhd_rules::magnetic_extension_rules(nd::uint depth, nd::uint block_size, rational::number_t iteration)
{
    return block_indexes(depth, block_size)
    | seq::map([=] (auto index) -> named_rule_t
    {
        auto bf = key(data_field::face_magnetic_flux_density).iteration(iteration);
        auto n1 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), block_extent(depth), mesh::axis_3d::i);
        auto n2 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), block_extent(depth), mesh::axis_3d::j);
        auto n3 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), block_extent(depth), mesh::axis_3d::k);

        auto neighbor_ids = seq::concat(n1, n2, n3)
        | seq::map([=] (auto i) { return multilevel_index_t{depth, mesh::to_numeric_array(i)}; })
        | seq::map(bf.bind_block())
        | seq::to<std::vector>();

        return rule(bf.extended().block(index), extend_face_magnetic_flux_density(block_size, 1), neighbor_ids);
    })
    | seq::to_dynamic();
}




//=============================================================================
seq::dynamic_sequence_t<named_rule_t> mhd_rules::godunov_data_rules(nd::uint depth, nd::uint block_size, rational::number_t iteration)
{
    return block_indexes(depth, block_size)
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
seq::dynamic_sequence_t<named_rule_t> mhd_rules::electromotive_force_rules(nd::uint depth, nd::uint block_size, rational::number_t iteration)
{
    return block_indexes(depth, block_size)
    | seq::map([=] (auto index) -> named_rule_t
    {
        auto ee = key(data_field::edge_electromotive_density).iteration(iteration).block(index);
        auto gf = key(data_field::face_godunov_data         ).iteration(iteration).block(index);
        auto fm = wrap<face_godunov_data_t>(electromotive_forces);
        return rule(ee, fm, {gf});
    })
    | seq::to_dynamic();
}




//=============================================================================
seq::dynamic_sequence_t<named_rule_t> mhd_rules::global_primitive_array_rules(nd::uint depth, nd::uint block_size, rational::number_t iteration)
{
    auto pc = key(data_field::cell_primitive_variables).iteration(iteration);
    auto arg_keys = block_indexes(depth, block_size) | seq::map(pc.bind_block()) | seq::to<std::vector>();

    auto tile = [=] (std::vector<product_t> block_vector) -> product_t
    {
        auto block_map = seq::view(block_vector)
        | seq::map([] (auto p) { return std::get<cell_primitive_variables_t>(p); })
        | seq::keys(nd::index_space(block_extent(depth)))
        | seq::to_dict<std::map>();

        return mesh::tile_blocks(block_map, block_extent(depth)) | nd::to_shared();
    };
    return seq::just(rule(pc, tile, arg_keys)) | seq::to_dynamic();
}
