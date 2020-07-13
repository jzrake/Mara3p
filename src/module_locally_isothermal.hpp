/**
 ==============================================================================
 Copyright 2020, Jonathan Zrake

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
#include <map>
#include "core_bqo_tree.hpp"
#include "core_bsp_tree.hpp"
#include "core_ndarray.hpp"
#include "mesh_geometry.hpp"
#include "parallel_computable.hpp"
#include "physics_iso2d.hpp"




//=============================================================================
namespace modules
{




//=============================================================================
namespace locally_isothermal
{

using namespace dimensional;
using coords_t              = mesh::cartesian_2d::coords_t;
using mesh_geometry_t       = mesh::cartesian_2d::geometry_t;
using mesh_topology_t       = bsp::shared_tree<mesh::block_index_t<2>, 4>;
using conserved_array_t     = nd::shared_array<iso2d::conserved_density_t, 2>;
using primitive_array_t     = nd::shared_array<iso2d::primitive_t, 2>;
using godunov_f_array_t     = nd::shared_array<iso2d::flux_vector_t, 2>;
using conserved_tree_t      = bsp::shared_tree<mpr::computable<conserved_array_t>, 4>;
using primitive_tree_t      = bsp::shared_tree<mpr::computable<primitive_array_t>, 4>;
using primitive_mapping_t   = std::function<iso2d::primitive_t(coords_t)>;
using temperature_mapping_t = std::function<dimensional::unit_specific_energy(coords_t)>;
using source_terms_t        = std::function<conserved_array_t(coords_t, unit_time, unit_time, conserved_array_t)>;
using source_terms_map_t    = std::map<std::string, source_terms_t>;
using source_totals_t       = std::map<std::string, iso2d::conserved_density_t>;
using unit_viscosity        = dimensional::quantity_t<0, 2,-1>;




//=============================================================================
conserved_array_t initial_conserved_array(
    mesh_geometry_t mesh_geometry,
    mesh::block_index_t<2> block,
    primitive_mapping_t initial);

conserved_tree_t initial_conserved_tree(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry,
    primitive_mapping_t initial);

primitive_array_t recover_primitive_array(
    conserved_array_t uc);

primitive_array_t estimate_gradient(
    primitive_array_t pc,
    unsigned long axis,
    double plm_theta);

conserved_array_t updated_conserved(
    conserved_array_t uc,
    primitive_array_t pe,
    unit_time time,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    mesh::block_index_t<2> block,
    unit_viscosity kinematic_viscosity,
    temperature_mapping_t temperature,
    double plm_theta);

conserved_tree_t updated_conserved_tree(
    conserved_tree_t conserved,
    unit_time time,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    temperature_mapping_t temperature,
    unit_viscosity kinematic_viscosity,
    double plm_theta);

unit_length smallest_cell_size(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry);

std::size_t total_zones(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry);

} // namespace locally_isothermal

} // namespace modules




//=============================================================================
namespace h5 {

class Group;

void read(const Group& group, std::string name, modules::locally_isothermal::conserved_tree_t& conserved);
void write(const Group& group, std::string name, const modules::locally_isothermal::conserved_tree_t& conserved);

} // namespace h5
