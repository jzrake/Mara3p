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




#include "app_serial_ndarray.hpp"
#include "app_serial_numeric_tuple.hpp"
#include "app_serial_std_tuple.hpp"
#include "core_memoize.hpp"
#include "core_ndarray_ops.hpp"
#include "mara.hpp"
#include "mesh_amr.hpp"
#include "module_locally_isothermal.hpp"
#include "scheme_plm_gradient.hpp"




using namespace std::placeholders;
using namespace modules;




#if MARA_COMPILE_LOCALLY_ISOTHERMAL // <---------------------------------------


//=============================================================================
locally_isothermal::conserved_array_t locally_isothermal::initial_conserved_array(
    mesh_geometry_t mesh_geometry,
    mesh::block_index_t<2> block,
    primitive_mapping_t initial)
{
    return mesh_geometry.cell_coordinates(block)
    | nd::map(initial)
    | nd::map(iso2d::conserved_density)
    | nd::to_shared();
}




//=============================================================================
locally_isothermal::conserved_tree_t locally_isothermal::initial_conserved_tree(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry,
    primitive_mapping_t initial)
{
    auto u0 = std::bind(initial_conserved_array, mesh_geometry, _1, initial);
    auto uc = [u0] (auto block) { return mpr::from(std::bind(u0, block)).name("U"); };
    return mesh_topology | bsp::maps(uc);
}




//=============================================================================
locally_isothermal::primitive_array_t locally_isothermal::recover_primitive_array(conserved_array_t uc)
{
    return uc
    | nd::map(iso2d::recover_primitive)
    | nd::to_shared();
}




//=============================================================================
locally_isothermal::primitive_array_t locally_isothermal::estimate_gradient(primitive_array_t pc, unsigned long axis, double plm_theta)
{
    return pc
    | nd::adjacent_zip3(axis)
    | nd::map(mara::plm_gradient(plm_theta))
    | nd::to_shared();
}




namespace modules::locally_isothermal {

//=============================================================================
template<typename PrimitiveArray, typename GradientArrayL, typename GradientArrayT>
godunov_f_array_t godunov_and_viscous_fluxes(
    mesh_geometry_t mesh_geometry,
    PrimitiveArray pc,
    GradientArrayL gc_long,
    GradientArrayT gc_tran,
    mesh::block_index_t<2> block,
    unit_viscosity kinematic_viscosity,
    temperature_mapping_t temperature,
    unsigned long axis)
{
    auto riemann = [axis, temperature] (auto pl, auto pr, auto xf)
    {
        auto cs2 = temperature(xf);
        return iso2d::riemann_hlle(pl, pr, cs2, geometric::unit_vector_on(axis + 1));
    };

    auto xf = mesh_geometry.face_coordinates(block, axis);
    auto pl = (pc - 0.5 * gc_long) | nd::select(axis, 1)     | nd::to_dynamic();
    auto pr = (pc + 0.5 * gc_long) | nd::select(axis, 0, -1) | nd::to_dynamic();
    auto fx = nd::zip(pr, pl, xf) | nd::mapv(riemann);

    auto viscous_flux = [=] ()
    {
        auto nu = kinematic_viscosity;
        auto sl = map(pl, iso2d::mass_density);
        auto sr = map(pr, iso2d::mass_density);
        auto mu = nu * 0.5 * (sl + sr);
        auto dl = mesh_geometry.cell_spacing(block);
        auto block_size = mesh_geometry.block_size;

        switch (axis)
        {
            case 0:
            {
                auto dx_vx = gc_long | nd::map(iso2d::velocity_1) | nd::adjacent_mean(0) | nd::divide(dl);
                auto dy_vx = gc_tran | nd::map(iso2d::velocity_1) | nd::adjacent_mean(0) | nd::divide(dl);
                auto dx_vy = gc_long | nd::map(iso2d::velocity_2) | nd::adjacent_mean(0) | nd::divide(dl);
                auto dy_vy = gc_tran | nd::map(iso2d::velocity_2) | nd::adjacent_mean(0) | nd::divide(dl);

                auto tauxs = nd::zeros<dimensional::quantity_t<1,-1,-1>>(block_size + 1, block_size);
                auto tauxx = mu * (dx_vx - dy_vy);
                auto tauxy = mu * (dx_vy + dy_vx);

                return nd::zip(tauxs, -tauxx, -tauxy) | nd::construct<iso2d::flux_vector_t>() | nd::to_dynamic();
            }
            case 1:
            {
                auto dx_vx = gc_tran | nd::map(iso2d::velocity_1) | nd::adjacent_mean(1) | nd::divide(dl);
                auto dy_vx = gc_long | nd::map(iso2d::velocity_1) | nd::adjacent_mean(1) | nd::divide(dl);
                auto dx_vy = gc_tran | nd::map(iso2d::velocity_2) | nd::adjacent_mean(1) | nd::divide(dl);
                auto dy_vy = gc_long | nd::map(iso2d::velocity_2) | nd::adjacent_mean(1) | nd::divide(dl);

                auto tauys = nd::zeros<dimensional::quantity_t<1,-1,-1>>(block_size, block_size + 1);
                auto tauyx = mu * (dx_vy + dy_vx);
                auto tauyy =-mu * (dx_vx - dy_vy);

                return nd::zip(tauys, -tauyx, -tauyy) | nd::construct<iso2d::flux_vector_t>() | nd::to_dynamic();
            }
        }
        throw;
    };

    if (kinematic_viscosity == unit_viscosity(0.0))
    {        
        return fx | nd::to_shared();
    }
    return (fx + viscous_flux()) | nd::to_shared();
}

} // modules::locally_isothermal 




//=============================================================================
locally_isothermal::conserved_array_t locally_isothermal::updated_conserved(
    conserved_array_t uc,
    primitive_array_t pe,
    unit_time time,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    mesh::block_index_t<2> block,
    unit_viscosity kinematic_viscosity,
    temperature_mapping_t temperature,
    double plm_theta)
{
    auto xc = mesh_geometry.cell_coordinates(block);
    auto dl = mesh_geometry.cell_spacing(block);
    auto gx = estimate_gradient(pe, 0, 2.0);
    auto gy = estimate_gradient(pe, 1, 2.0);
    auto pc = pe | nd::select(0, 2, -2) | nd::select(1, 2, -2) | nd::to_shared();
    auto pe_x = pe | nd::select(1, 2, -2) | nd::select(0, 1, -1);
    auto gx_x = gx | nd::select(1, 2, -2);
    auto gy_x = gy | nd::select(1, 1, -1) | nd::select(0, 1, -1);
    auto pe_y = pe | nd::select(0, 2, -2) | nd::select(1, 1, -1);
    auto gy_y = gy | nd::select(0, 2, -2);
    auto gx_y = gx | nd::select(0, 1, -1) | nd::select(1, 1, -1);

    auto fx = godunov_and_viscous_fluxes(mesh_geometry, pe_x, gx_x, gy_x, block, kinematic_viscosity, temperature, 0);
    auto fy = godunov_and_viscous_fluxes(mesh_geometry, pe_y, gy_y, gx_y, block, kinematic_viscosity, temperature, 1);
    auto dfx = fx | nd::adjacent_diff(0);
    auto dfy = fy | nd::adjacent_diff(1);
    auto u1 = uc - (dfx + dfy) * dt / dl;

    return u1 | nd::to_shared();
}




//=============================================================================
locally_isothermal::conserved_tree_t locally_isothermal::updated_conserved_tree(
    conserved_tree_t conserved,
    unit_time time,
    unit_time dt,
    mesh_geometry_t mesh_geometry,
    temperature_mapping_t temperature,
    unit_viscosity kinematic_viscosity,
    double plm_theta)
{
    auto F = std::bind(updated_conserved, _1, _2, time, dt, mesh_geometry, _3, kinematic_viscosity, temperature, plm_theta);
    auto U = [F] (auto b) { return std::bind(F, _1, _2, b); };

    auto mesh  = indexes(conserved);
    auto uc    = conserved;
    auto pc    = uc   | bsp::maps(mpr::map(recover_primitive_array, "P"));
    auto pc_at = memoize_not_thread_safe<mesh::block_index_t<2>>([pc] (auto block) { return amr::get_or_create_block(pc, block).name("Pc-g"); });
    auto pe    = mesh | bsp::maps([pc, pc_at] (auto b) { return amr::extend_block(pc_at, b).name("Pe"); });
    auto u1    = mesh | bsp::maps([uc, pe, U] (auto b) { return zip(value_at(uc, b), value_at(pe, b)) | mpr::mapv(U(b), "U"); });

    return u1;
}




//=============================================================================
dimensional::unit_length locally_isothermal::smallest_cell_size(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry)
{
    return reduce(mesh_topology, [mesh_geometry] (auto seed, auto block)
    {
        return std::min(seed, mesh_geometry.cell_spacing(block));
    }, unit_length(std::numeric_limits<double>::max()));
}




//=============================================================================
std::size_t locally_isothermal::total_zones(
    mesh_topology_t mesh_topology,
    mesh_geometry_t mesh_geometry)
{
    return size(mesh_topology) * mesh_geometry.cells_per_block();
}




//=============================================================================
locally_isothermal::solution_t locally_isothermal::weighted_sum(solution_t s, solution_t t, rational::number_t b)
{
    return {
        s.iteration  *         b  + t.iteration *       (1 - b),
        s.time       *  double(b) + t.time      * double(1 - b),
        amr::weighted_sum_tree(s.conserved, t.conserved, b),
    };
}

#endif // MARA_COMPILE_LOCALLY_ISOTHERMAL
