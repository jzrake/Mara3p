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
#include "core_bsp_tree.hpp"
#include "core_bqo_tree.hpp"
#include "core_ndarray.hpp"
#include "physics_euler.hpp"
#include "parallel_computable.hpp"




//=============================================================================
namespace modules
{

namespace euler2d
{




using namespace dimensional;
using coords_t            = numeric::array_t<unit_length, 2>;
using mesh_topology_t     = bsp::shared_tree<bsp::tree_index_t<2>, 4>;
using conserved_array_t   = nd::shared_array<euler::conserved_density_t, 2>;
using primitive_array_t   = nd::shared_array<euler::primitive_t, 2>;
using conserved_tree_t    = bsp::shared_tree<mpr::computable<conserved_array_t>, 4>;
using primitive_tree_t    = bsp::shared_tree<mpr::computable<primitive_array_t>, 4>;
using primitive_mapping_t = std::function<euler::primitive_t(coords_t)>;




//=============================================================================
struct side_effect_t
{
    unit_time next_due = 0.0;
    int count = 0;
};

struct schedule_t
{
    side_effect_t checkpoint;
    side_effect_t time_series;
};

struct solution_t
{
    rational::number_t iteration;
    unit_time          time;
    conserved_tree_t   conserved;
};




//=============================================================================
class mesh_geometry_t
{
public:
    mesh_geometry_t() {}
    mesh_geometry_t(unit_length domain_size, int block_size);

    nd::shared_array<coords_t, 2> vert_coordinates(bsp::tree_index_t<2> block) const;
    nd::shared_array<coords_t, 2> face_coordinates(bsp::tree_index_t<2> block, bsp::uint axis) const;
    nd::shared_array<coords_t, 2> cell_coordinates(bsp::tree_index_t<2> block) const;
    unit_length cell_spacing(bsp::tree_index_t<2> block) const;
    std::size_t cells_per_block() const;

private:
    unit_length domain_size = 1.0;
    int block_size = 64;
};




//=============================================================================
conserved_array_t initial_conserved_array(
    mesh_geometry_t mesh_geometry,
    bsp::tree_index_t<2> block,
    primitive_mapping_t initial,
    double gamma_law_index);

conserved_tree_t initial_conserved_tree(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry,
    primitive_mapping_t initial,
    double gamma_law_index);

primitive_array_t recover_primitive_array(
    conserved_array_t uc,
    double gamma_law_index);

primitive_array_t estimate_gradient(
    primitive_array_t pc,
    bsp::uint axis,
    double plm_theta);

conserved_array_t updated_conserved(
    conserved_array_t uc,
    primitive_array_t pe,
    unit_time time,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    bsp::tree_index_t<2> block,
    double plm_theta,
    double gamma_law_index);

solution_t updated_solution(
    solution_t solution,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    double plm_theta,
    double gamma_law_index);

unit_length smallest_cell_size(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry);

std::size_t total_cells(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry);

solution_t weighted_sum(solution_t s, solution_t t, rational::number_t b);



} // namespace euler2d

} // namespace modules




//=============================================================================
namespace h5 {

class Group;

void read(const Group& group, std::string name, modules::euler2d::conserved_tree_t& conserved);
void read(const Group& group, std::string name, modules::euler2d::side_effect_t& side_effect);
void read(const Group& group, std::string name, modules::euler2d::schedule_t& schedule);
void read(const Group& group, std::string name, modules::euler2d::solution_t& solution);

void write(const Group& group, std::string name, const modules::euler2d::conserved_tree_t& conserved);
void write(const Group& group, std::string name, const modules::euler2d::side_effect_t& side_effect);
void write(const Group& group, std::string name, const modules::euler2d::schedule_t& schedule);
void write(const Group& group, std::string name, const modules::euler2d::solution_t& solution);

} // namespace h5
