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
namespace schemes
{

namespace euler2d
{




using namespace dimensional;
using coords_t            = numeric::array_t<unit_length, 2>;
using conserved_array_t   = nd::shared_array<euler::conserved_density_t, 2>;
using primitive_array_t   = nd::shared_array<euler::primitive_t, 2>;
using godunov_f_array_t   = nd::shared_array<euler::flux_vector_t, 2>;
using conserved_tree_t    = bsp::shared_tree<mpr::computable<conserved_array_t>, 4>;
using primitive_mapping_t = std::function<euler::primitive_t(coords_t)>;
using unit_viscosity      = dimensional::quantity_t<0, 2,-1>;




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
class mesh_coordinates_t
{
public:
    mesh_coordinates_t() {}
    mesh_coordinates_t(unit_length domain_size, int block_size);

    nd::shared_array<coords_t, 2> vert_coordinates(bsp::tree_index_t<2> block) const;
    nd::shared_array<coords_t, 2> face_coordinates(bsp::tree_index_t<2> block, bsp::uint axis) const;
    nd::shared_array<coords_t, 2> cell_coordinates(bsp::tree_index_t<2> block) const;
    unit_length cell_size(bsp::tree_index_t<2> block) const;

private:
    unit_length domain_size = 1.0;
    int block_size = 64;
};




//=============================================================================
class scheme_t
{
public:
    conserved_array_t initial_conserved_array(bsp::tree_index_t<2> block, primitive_mapping_t initial) const;
    primitive_array_t recover_primitive_array(conserved_array_t uc) const;
    primitive_array_t estimate_gradient(primitive_array_t pc, bsp::uint axis) const;
    conserved_array_t updated_conserved(
        conserved_array_t uc,
        primitive_array_t pe,
        unit_time time,
        unit_time dt,
        bsp::tree_index_t<2> block) const;
private:
    double plm_theta = 2.0;
    double gamma_law_index = 5. / 3;
    mesh_coordinates_t mesh;
};




} // namespace euler2d

} // namespace schemes
