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




#pragma once
#include <variant>
#include "core_rational.hpp"
#include "core_sequence.hpp"
#include "parallel_dependency_graph.hpp"
#include "physics_mhd.hpp"
#include "scheme_mhd_v2.hpp"




//=============================================================================
namespace mhd_rules
{




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
using product_t = std::variant<
    mhd_scheme_v2::cell_primitive_variables_t,
    mhd_scheme_v2::cell_conserved_density_t,
    mhd_scheme_v2::face_godunov_data_t,
    mhd_scheme_v2::face_magnetic_flux_density_t,
    mhd_scheme_v2::edge_electromotive_density_t>;

using product_identifier_t = std::tuple<
    rational::number_t,
    mhd_scheme_v2::multilevel_index_t,
    data_field,
    extended_status>;

using DependencyGraph = mara::DependencyGraph<product_identifier_t, product_t>;
using mapping_type    = DependencyGraph::mapping_type;




//=============================================================================
using named_rule_t = std::tuple<product_identifier_t, mapping_type, std::vector<product_identifier_t>>;

inline named_rule_t rule(product_identifier_t key, mapping_type mapping, std::vector<product_identifier_t> args)
{
    return named_rule_t{key, mapping, args};
}




//=============================================================================
seq::dynamic_sequence_t<named_rule_t> initial_condition_rules     (nd::uint depth, nd::uint block_size);
seq::dynamic_sequence_t<named_rule_t> recover_primitive_rules     (nd::uint depth, nd::uint block_size, rational::number_t iteration);
seq::dynamic_sequence_t<named_rule_t> primitive_extension_rules   (nd::uint depth, nd::uint block_size, rational::number_t iteration);
seq::dynamic_sequence_t<named_rule_t> magnetic_extension_rules    (nd::uint depth, nd::uint block_size, rational::number_t iteration);
seq::dynamic_sequence_t<named_rule_t> godunov_data_rules          (nd::uint depth, nd::uint block_size, rational::number_t iteration);
seq::dynamic_sequence_t<named_rule_t> electromotive_force_rules   (nd::uint depth, nd::uint block_size, rational::number_t iteration);
seq::dynamic_sequence_t<named_rule_t> global_primitive_array_rules(nd::uint depth, nd::uint block_size, rational::number_t iteration);




//=============================================================================
struct key_factory_t
{


    //=============================================================================
    key_factory_t iteration(rational::number_t v)                        const { auto n = id; std::get<0>(n) = v; return {n}; }
    key_factory_t block    (mhd_scheme_v2::multilevel_index_t v)         const { auto n = id; std::get<1>(n) = v; return {n}; }
    key_factory_t field    (data_field v)                                const { auto n = id; std::get<2>(n) = v; return {n}; }
    key_factory_t extended (extended_status v=extended_status::extended) const { auto n = id; std::get<3>(n) = v; return {n}; }


    /**
     * @brief      Product_identifier_t user-conversion operator.
     */
    operator product_identifier_t() const
    {
        return id;
    }


    /**
     * @brief      Return a function that returns a product identifier with the
     *             same values as this factory object, but with the block index
     *             given by the argument to the function.
     *
     * @return     A function (multilevel_index_t) -> product_identifier_t
     */
    auto bind_block() const
    {
        return [id=id] (mhd_scheme_v2::multilevel_index_t v)
        {
            return key_factory_t{id}.block(v).id;
        };
    }

    product_identifier_t id;
};




//=============================================================================
inline key_factory_t key(data_field field)
{
    return key_factory_t().field(field);
}

inline key_factory_t key(rational::number_t iteration)
{
    return key_factory_t().iteration(iteration);
}

} // namespace mhd_rules
