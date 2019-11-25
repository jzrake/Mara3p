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




#include "scheme_mhd.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"
#include "physics_mhd.hpp"
#include "mesh_cartesian_3d.hpp"




//=============================================================================
using namespace util;
using namespace std::placeholders;
using geometric::unit_vector_on;
static const double gamma_law_index = 5. / 3;




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
template<typename P>
static auto edge_emf_from_face(nd::array_t<P, 3> ef1, nd::array_t<P, 3> ef2, nd::array_t<P, 3> ef3)
{
    auto c1 = nd::map(geometric::component(1));
    auto c2 = nd::map(geometric::component(2));
    auto c3 = nd::map(geometric::component(3));
    auto a1 = nd::adjacent_mean(0);
    auto a2 = nd::adjacent_mean(1);
    auto a3 = nd::adjacent_mean(2);

    auto ee1a = ef2 | c1 | a3;
    auto ee2a = ef3 | c2 | a1;
    auto ee3a = ef1 | c3 | a2;
    auto ee1b = ef3 | c1 | a2;
    auto ee2b = ef1 | c2 | a3;
    auto ee3b = ef2 | c3 | a1;

    auto ee1 = 0.5 * (ee1a + ee1b);
    auto ee2 = 0.5 * (ee2a + ee2b);
    auto ee3 = 0.5 * (ee3a + ee3b);

    return std::tuple(ee1, ee2, ee3);
}




template<typename P, typename Q>
static auto godunov_fluxes(nd::array_t<P, 3> pc, nd::array_t<Q, 3> bf1, nd::array_t<Q, 3> bf2, nd::array_t<Q, 3> bf3)
{
    auto pf1 = nd::zip(pc | nd::select(0, 0, -1), pc | nd::select(0, 1), bf1);
    auto pf2 = nd::zip(pc | nd::select(1, 0, -1), pc | nd::select(1, 1), bf2);
    auto pf3 = nd::zip(pc | nd::select(2, 0, -1), pc | nd::select(2, 1), bf3);

    auto [ff1, ef1] = nd::unzip(pf1 | nd::map(rs1) | nd::to_shared());
    auto [ff2, ef2] = nd::unzip(pf2 | nd::map(rs2) | nd::to_shared());
    auto [ff3, ef3] = nd::unzip(pf3 | nd::map(rs3) | nd::to_shared());

    return std::tuple(ff1, ff2, ff3, ef1, ef2, ef3);
}




//=============================================================================
nd::shared_array<mhd::primitive_t, 3> primitive_array(
    nd::shared_array<mhd::conserved_density_t, 3> uc,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf1,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf2,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf3,
    double gamma_law_index)
{
    auto bc = mesh::face_to_cell(bf1, bf2, bf3);
    return zip(uc, bc) | nd::map(apply_to(std::bind(mhd::recover_primitive, _1, _2, gamma_law_index))) | nd::to_shared();
}




//=============================================================================
std::tuple<
    nd::shared_array<mhd::flux_vector_t, 3>,
    nd::shared_array<mhd::flux_vector_t, 3>,
    nd::shared_array<mhd::flux_vector_t, 3>,
    nd::shared_array<mhd::unit_electric_field, 3>,
    nd::shared_array<mhd::unit_electric_field, 3>,
    nd::shared_array<mhd::unit_electric_field, 3>>
flux_arrays(
    nd::shared_array<mhd::primitive_t, 3> pc,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf1,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf2,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf3)
{
    auto [ff1, ff2, ff3, ef1, ef2, ef3] = godunov_fluxes(pc, bf1, bf2, bf3);
    auto [ee1, ee2, ee3] = edge_emf_from_face(ef1, ef2, ef3);
    auto ev = nd::to_shared();
    return std::tuple(ff1|rt1|ev, ff2|rt2|ev, ff3|rt3|ev, ee1|rl1|ev, ee2|rl2|ev, ee3|rl3|ev);
}




//=============================================================================
std::tuple<
    dimensional::unit_time,
    nd::shared_array<mhd::conserved_density_t, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>>
advance(dimensional::unit_time time,
        nd::shared_array<mhd::conserved_density_t, 3> uc,
        nd::shared_array<mhd::unit_magnetic_field, 3> bf1,
        nd::shared_array<mhd::unit_magnetic_field, 3> bf2,
        nd::shared_array<mhd::unit_magnetic_field, 3> bf3,
        dimensional::unit_length dl)
{
    auto pc = primitive_array(uc, bf1, bf2, bf3, gamma_law_index);

    auto pce  = pc  | mesh::extend_periodic(1) | nd::to_shared();
    auto bf1e = bf1 | nd::extend_periodic(1) | nd::extend_periodic(2) | nd::to_shared();
    auto bf2e = bf2 | nd::extend_periodic(2) | nd::extend_periodic(0) | nd::to_shared();
    auto bf3e = bf3 | nd::extend_periodic(0) | nd::extend_periodic(1) | nd::to_shared();

    auto [ff1, ff2, ff3, ee1, ee2, ee3] = flux_arrays(pce, bf1e, bf2e, bf3e);
    auto [cf1, cf2, cf3] = mesh::solenoidal_difference(ee1, ee2, ee3);
    auto dc              = mesh::divergence_difference(ff1, ff2, ff3);

    auto dt   = 0.33 * dl / max(pc | nd::map(std::bind(mhd::fastest_speed, _1, gamma_law_index)));
    auto uc_  = uc  - dc  * dt / dl;
    auto bf1_ = bf1 - cf1 * dt / dl;
    auto bf2_ = bf2 - cf2 * dt / dl;
    auto bf3_ = bf3 - cf3 * dt / dl;

    return std::tuple(
        time + dt,
        uc_  | nd::to_shared(),
        bf1_ | nd::to_shared(),
        bf2_ | nd::to_shared(),
        bf3_ | nd::to_shared());
}
