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




#include "core_memoize.hpp"
#include "core_ndarray_ops.hpp"
#include "module_euler.hpp"
#include "mesh_amr.hpp"
#include "scheme_plm_gradient.hpp"




using namespace std::placeholders;
using namespace modules;
using namespace euler2d;




//=============================================================================
conserved_array_t euler2d::initial_conserved_array(
    mesh_geometry_t mesh_geometry,
    mesh::block_index_t<2> block,
    primitive_mapping_t initial,
    double gamma_law_index)
{
    return mesh_geometry.cell_coordinates(block)
    | nd::map(initial)
    | nd::map(std::bind(euler::conserved_density, _1, gamma_law_index))
    | nd::to_shared();
}




//=============================================================================
conserved_tree_t euler2d::initial_conserved_tree(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry,
    primitive_mapping_t initial,
    double gamma_law_index)
{
    auto u0 = std::bind(initial_conserved_array, mesh_geometry, _1, initial, gamma_law_index);
    auto uc = [u0] (auto block) { return mpr::from(std::bind(u0, block)).name("U"); };
    return mesh_topology | bsp::maps(uc);
}




//=============================================================================
primitive_array_t euler2d::recover_primitive_array(conserved_array_t uc, double gamma_law_index)
{
    return uc
    | nd::map(std::bind(euler::recover_primitive, _1, gamma_law_index))
    | nd::to_shared();
}




//=============================================================================
primitive_array_t euler2d::estimate_gradient(primitive_array_t pc, unsigned long axis, double plm_theta)
{
    return pc
    | nd::adjacent_zip3(axis)
    | nd::map(mara::plm_gradient(plm_theta))
    | nd::to_shared();
}




//=============================================================================
conserved_array_t euler2d::updated_conserved(
    conserved_array_t uc,
    primitive_array_t pe,
    unit_time time,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    mesh::block_index_t<2> block,
    double plm_theta,
    double gamma_law_index)
{
    auto dl = mesh_geometry.cell_spacing(block);
    auto gx = estimate_gradient(pe, 0, plm_theta);
    auto gy = estimate_gradient(pe, 1, plm_theta);
    auto pc = pe | nd::select(0, 2, -2) | nd::select(1, 2, -2) | nd::to_shared();

    auto pe_x = pe | nd::select(1, 2, -2) | nd::select(0, 1, -1);
    auto gx_x = gx | nd::select(1, 2, -2);
    auto pe_y = pe | nd::select(0, 2, -2) | nd::select(1, 1, -1);
    auto gy_y = gy | nd::select(0, 2, -2);

    auto riemann_x = std::bind(euler::riemann_hlle, _1, _2, geometric::unit_vector_on(1), gamma_law_index);
    auto riemann_y = std::bind(euler::riemann_hlle, _1, _2, geometric::unit_vector_on(2), gamma_law_index);

    auto pl_x = (pe_x - 0.5 * gx_x) | nd::select(0, 1);
    auto pr_x = (pe_x + 0.5 * gx_x) | nd::select(0, 0, -1);
    auto pl_y = (pe_y - 0.5 * gy_y) | nd::select(1, 1);
    auto pr_y = (pe_y + 0.5 * gy_y) | nd::select(1, 0, -1);

    auto fx = nd::zip(pr_x, pl_x) | nd::mapv(riemann_x);
    auto fy = nd::zip(pr_y, pl_y) | nd::mapv(riemann_y);

    auto dfx = fx | nd::adjacent_diff(0);
    auto dfy = fy | nd::adjacent_diff(1);

    return (uc - (dfx + dfy) * dt / dl) | nd::to_shared();
}




//=============================================================================
solution_t euler2d::updated_solution(
    solution_t solution,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    double plm_theta,
    double gamma_law_index)
{


    auto F = std::bind(updated_conserved, _1, _2, solution.time, dt, mesh_geometry, _3, plm_theta, gamma_law_index);
    auto U = [F] (auto b) { return std::bind(F, _1, _2, b); };

    auto mesh  = indexes(solution.conserved);
    auto uc    = solution.conserved;
    auto pc    = uc   | bsp::maps(mpr::map([g=gamma_law_index] (auto u) { return recover_primitive_array(u, g); }, "P"));
    auto pc_at = memoize_not_thread_safe<mesh::block_index_t<2>>([pc] (auto block) { return amr::get_or_create_block(pc, block); });
    auto pe    = mesh | bsp::maps([pc, pc_at] (auto b) { return amr::extend_block(pc_at, b); });
    auto u1    = mesh | bsp::maps([uc, pe, U] (auto b) { return zip(value_at(uc, b), value_at(pe, b)) | mpr::mapv(U(b)); });

    return solution_t{solution.iteration + 1, solution.time + dt, u1};
}




//=============================================================================
unit_length euler2d::smallest_cell_size(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry)
{
    return reduce(mesh_topology, [mesh_geometry] (auto seed, auto block)
    {
        return std::min(seed, mesh_geometry.cell_spacing(block));
    }, unit_length(std::numeric_limits<double>::max()));
}




//=============================================================================
std::size_t euler2d::total_cells(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry)
{
    return size(mesh_topology) * mesh_geometry.cells_per_block();
}




//=============================================================================
solution_t euler2d::weighted_sum(solution_t s, solution_t t, rational::number_t b)
{
    return {
        s.iteration  *         b  + t.iteration *       (1 - b),
        s.time       *  double(b) + t.time      * double(1 - b),
        amr::weighted_sum_tree(s.conserved, t.conserved, b),
    };
}
