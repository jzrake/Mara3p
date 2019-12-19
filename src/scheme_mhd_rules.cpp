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




#include <variant>
#include "core_sequence.hpp"
#include "scheme_mhd_rules.hpp"




using namespace mhd_scheme_v2;
using namespace mhd_rules;




//=============================================================================
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
std::array<named_rule_t, 3> mhd_rules::initial_condition(multilevel_index_t index, nd::uint block_size, mhd_scheme_v2::primitive_function_t p, mhd_scheme_v2::vector_potential_function_t a)
{
    auto nb = unsigned(1 << index.level);
    auto dl = dimensional::unit_length(1.0 / block_size / nb);

    auto k_ae = key(data_field::edge_vector_potential)     .block(index);
    auto k_bf = key(data_field::face_magnetic_flux_density).block(index);
    auto k_uc = key(data_field::cell_conserved_density)    .block(index);

    auto f_ae = wrap([=] () { return construct_vector_potential(index, block_size, a); });
    auto f_bf = wrap<edge_electromotive_density_t>([=] (auto ee) { return curl(ee, dl); });
    auto f_uc = wrap<face_magnetic_flux_density_t>([=] (auto bf) { return construct_conserved(index, block_size, bf, p); });

    auto r_ae = rule(k_ae, f_ae, {});
    auto r_bf = rule(k_bf, f_bf, {k_ae});
    auto r_uc = rule(k_uc, f_uc, {k_bf});

    return std::array{r_ae, r_bf, r_uc};
}




//=============================================================================
named_rule_t mhd_rules::recover_primitive(multilevel_index_t index, rational::number_t iteration)
{
    auto bf = key(iteration).block(index).field(data_field::face_magnetic_flux_density);
    auto uc = key(iteration).block(index).field(data_field::cell_conserved_density);
    auto pc = key(iteration).block(index).field(data_field::cell_primitive_variables);
    auto fm = wrap<cell_conserved_density_t, face_magnetic_flux_density_t>(mhd_scheme_v2::primitive_array);
    return rule(pc, fm, {uc, bf});
}




//=============================================================================
named_rule_t mhd_rules::extend_primitive(multilevel_index_t index, rational::number_t iteration, nd::uint block_size)
{
    auto pc = key(data_field::cell_primitive_variables).iteration(iteration);

    auto neighbor_ids = mesh::neighbors_27(mesh::to_uivec(index.coordinates), mesh::block_extent(index.level))
    | seq::map([=] (auto i) { return multilevel_index_t{index.level, mesh::to_numeric_array(i)}; })
    | seq::map(pc.bind_block())
    | seq::to<std::vector>();

    return rule(pc.extended().block(index), extend_cell_primitive_variables(block_size, 1), neighbor_ids);
}




//=============================================================================
named_rule_t mhd_rules::extend_magnetic(multilevel_index_t index, rational::number_t iteration, nd::uint block_size)
{
    auto bf = key(data_field::face_magnetic_flux_density).iteration(iteration);
    auto n1 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), mesh::block_extent(index.level), mesh::axis_3d::i);
    auto n2 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), mesh::block_extent(index.level), mesh::axis_3d::j);
    auto n3 = mesh::neighbors_9(mesh::to_uivec(index.coordinates), mesh::block_extent(index.level), mesh::axis_3d::k);

    auto neighbor_ids = seq::concat(n1, n2, n3)
    | seq::map([=] (auto i) { return multilevel_index_t{index.level, mesh::to_numeric_array(i)}; })
    | seq::map(bf.bind_block())
    | seq::to<std::vector>();

    return rule(bf.extended().block(index), extend_face_magnetic_flux_density(block_size, 1), neighbor_ids);
}




//=============================================================================
named_rule_t mhd_rules::godunov_data(multilevel_index_t index, rational::number_t iteration)
{
    auto pc = key(data_field::cell_primitive_variables)  .iteration(iteration).block(index).extended();
    auto bf = key(data_field::face_magnetic_flux_density).iteration(iteration).block(index).extended();
    auto gf = key(data_field::face_godunov_data)         .iteration(iteration).block(index);
    auto fm = wrap<cell_primitive_variables_t, face_magnetic_flux_density_t>(godunov_fluxes);
    return rule(gf, fm, {pc, bf});
}




//=============================================================================
named_rule_t mhd_rules::electromotive_force(multilevel_index_t index, rational::number_t iteration)
{
    auto ee = key(data_field::edge_electromotive_density).iteration(iteration).block(index);
    auto gf = key(data_field::face_godunov_data         ).iteration(iteration).block(index);
    auto fm = wrap<face_godunov_data_t>(electromotive_forces);
    return rule(ee, fm, {gf});
}




//=============================================================================
named_rule_t mhd_rules::update_conserved(multilevel_index_t index, rational::number_t iteration, dimensional::unit_time dt)
{
    auto gf0 = key(data_field::face_godunov_data)     .block(index).iteration(iteration - 1);
    auto uc0 = key(data_field::cell_conserved_density).block(index).iteration(iteration - 1);
    auto uc1 = key(data_field::cell_conserved_density).block(index).iteration(iteration);

    auto fm = [dt, index] (std::vector<product_t> args)
    {
        auto uc = std::get<cell_conserved_density_t>(args.at(0));
        auto gf = std::get<face_godunov_data_t>     (args.at(1));
        auto dl = dimensional::unit_length(1.0) / double(1 << index.level);

        return updated_conserved_density(uc, gf, dt, dl);
    };
    return rule(uc1, fm, {uc0, gf0});
}




//=============================================================================
named_rule_t mhd_rules::update_magnetic(multilevel_index_t index, rational::number_t iteration, dimensional::unit_time dt)
{
    auto ee0 = key(data_field::edge_electromotive_density).block(index).iteration(iteration - 1);
    auto bf0 = key(data_field::face_magnetic_flux_density).block(index).iteration(iteration - 1);
    auto bf1 = key(data_field::face_magnetic_flux_density).block(index).iteration(iteration);

    auto fm = [dt, index] (std::vector<product_t> args)
    {
        auto bf = std::get<face_magnetic_flux_density_t>(args.at(0));
        auto ee = std::get<edge_electromotive_density_t>(args.at(1));
        auto dl = dimensional::unit_length(1.0) / double(1 << index.level);

        return updated_magnetic_flux_density(bf, ee, dt, dl);
    };
    return rule(bf1, fm, {bf0, ee0});
}
