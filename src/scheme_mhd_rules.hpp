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
#include <string>
#include <variant>
#include "core_rational.hpp"
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

using multilevel_index_t   = mhd_scheme_v2::multilevel_index_t;
using product_mapping_t    = std::function<product_t(std::vector<product_t>)>;
using product_identifier_t = std::tuple<rational::number_t, multilevel_index_t, data_field, extended_status>;
using named_rule_t         = std::tuple<product_identifier_t, product_mapping_t, std::vector<product_identifier_t>>;




//=============================================================================
inline named_rule_t rule(product_identifier_t key, product_mapping_t mapping, std::vector<product_identifier_t> args)
{
    return named_rule_t{key, mapping, args};
}




//=============================================================================
std::array<named_rule_t, 3> initial_condition(multilevel_index_t index, nd::uint block_size, mhd_scheme_v2::primitive_function_t, mhd_scheme_v2::vector_potential_function_t);
named_rule_t recover_primitive  (multilevel_index_t index, rational::number_t iteration);
named_rule_t extend_primitive   (multilevel_index_t index, rational::number_t iteration, nd::uint block_size);
named_rule_t extend_magnetic    (multilevel_index_t index, rational::number_t iteration, nd::uint block_size);
named_rule_t godunov_data       (multilevel_index_t index, rational::number_t iteration);
named_rule_t electromotive_force(multilevel_index_t index, rational::number_t iteration);
named_rule_t update_conserved   (multilevel_index_t index, rational::number_t iteration, dimensional::unit_time dt);
named_rule_t update_magnetic    (multilevel_index_t index, rational::number_t iteration, dimensional::unit_time dt);




//=============================================================================
struct key_factory_t
{


    //=============================================================================
    key_factory_t iteration(rational::number_t v)                        const { auto n = id; std::get<0>(n) = v; return {n}; }
    key_factory_t block    (multilevel_index_t v)                        const { auto n = id; std::get<1>(n) = v; return {n}; }
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
        return [id=id] (multilevel_index_t v)
        {
            return key_factory_t{id}.block(v).id;
        };
    }


    //=============================================================================
    product_identifier_t id;
};




//=============================================================================
inline bool operator==(key_factory_t factory, product_identifier_t id)
{
    return factory.id == id;
}

inline bool operator==(product_identifier_t id, key_factory_t factory)
{
    return id == factory.id;
}

inline key_factory_t key(data_field field)
{
    return key_factory_t().field(field);
}

inline key_factory_t key(rational::number_t iteration)
{
    return key_factory_t().iteration(iteration);
}




//=============================================================================
inline std::string to_string(multilevel_index_t i)
{
    return std::to_string(i.level)
    + ":["
    + std::to_string(i.coordinates[0]) + ", "
    + std::to_string(i.coordinates[1]) + ", "
    + std::to_string(i.coordinates[2]) + "]";
}

inline std::string to_string(data_field f)
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

inline std::string to_string(extended_status s)
{
    switch (s)
    {
        case extended_status::not_extended:  return " [] ";
        case extended_status::extended:      return "[[]]";
        default: return "";
    }
}

inline std::string to_string(product_identifier_t id)
{
    auto [a, b, c, d] = id;
    return
      to_string(a) + " - "
    + to_string(b) + " - "
    + to_string(c) + " - "
    + to_string(d);
}

} // namespace mhd_rules
