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




#include "scheme_mhd_v2.hpp"




//=============================================================================
using namespace std::placeholders;
using namespace mhd_scheme_v2;
using util::apply_to;
using geometric::unit_vector_on;




static auto gamma_law_index = 5. / 3;




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
cell_primitive_variables_t mhd_scheme_v2::primitive_array(cell_conserved_density_t uc, face_magnetic_flux_density_t bf)
{
    auto c2p = apply_to(std::bind(mhd::recover_primitive, _1, _2, gamma_law_index));
    auto [bf1, bf2, bf3] = bf;
    auto bc = mesh::face_to_cell(bf1, bf2, bf3);
    return zip(uc, bc) | nd::map(c2p) | nd::to_shared();
}




//=============================================================================
edge_electromotive_density_t mhd_scheme_v2::construct_vector_potential(
    multilevel_index_t index,
    unsigned block_size,
    vector_potential_function_t vector_potential)
{
    auto A = [vector_potential] (unsigned dir) { return util::compose(geometric::component(dir), vector_potential); };
    auto dt = dimensional::unit_time(1.0);
    auto [xe1, xe2, xe3] = mesh::edge_positions(mesh::construct_vertices<dimensional::unit_length>(index.level, index.coordinates, block_size));

    return {
        ((xe1 | nd::map(A(1))) / dt) | nd::to_shared(),
        ((xe2 | nd::map(A(2))) / dt) | nd::to_shared(),
        ((xe3 | nd::map(A(3))) / dt) | nd::to_shared(),
    };
}




//=============================================================================
cell_conserved_density_t mhd_scheme_v2::construct_conserved(
    multilevel_index_t index,
    unsigned block_size,
    face_magnetic_flux_density_t bf,
    primitive_function_t primitive_function)
{
    auto p2c = std::bind(mhd::conserved_density, _1, gamma_law_index);
    auto [bf1, bf2, bf3] = bf;
    auto xc = mesh::cell_positions(mesh::construct_vertices<dimensional::unit_length>(index.level, index.coordinates, block_size));
    auto bc = mesh::face_to_cell(bf1, bf2, bf3);
    auto pc = nd::zip(xc, bc) | nd::map(apply_to(primitive_function));
    auto uc = pc | nd::map(p2c);
    return uc | nd::to_shared();
};




//=============================================================================
cell_primitive_variables_t mhd_scheme_v2::extend_cell_primitive_variables(std::vector<cell_primitive_variables_t> uc_vector, nd::uint block_size, nd::uint count)
{
    auto items = seq::view(uc_vector)
    | seq::keys(nd::index_space(3, 3, 3))
    | seq::to_dict<std::map>();

    return mesh::tile_blocks(items, nd::uivec(3, 3, 3))
    | mesh::remove_surface(block_size - count)
    | nd::to_shared();
}




//=============================================================================
face_magnetic_flux_density_t mhd_scheme_v2::extend_face_magnetic_flux_density(std::vector<face_magnetic_flux_density_t> bf_vector, nd::uint block_size, nd::uint count)
{
    auto bfs = seq::view(bf_vector) | seq::chunk(9);
    auto axs = seq::from(mesh::axis_3d::i, mesh::axis_3d::j, mesh::axis_3d::k);
    auto shp = seq::from(nd::uivec(1, 3, 3), nd::uivec(3, 1, 3), nd::uivec(3, 3, 1));

    auto extend = [r=block_size - count] (unsigned a, auto bf, mesh::axis_3d axis, nd::uivec_t<3> shape)
    {
        auto items = seq::adapt(bf)
        | seq::map([a] (auto x) { return x.at(a); })
        | seq::keys(nd::index_space(shape))
        | seq::to_dict<std::map>();

        return mesh::tile_blocks(items, shape)
        | mesh::remove_transverse(r, axis)
        | nd::to_shared();
    };

    return seq::zip(seq::range(3), bfs, axs, shp)
    | seq::map(util::apply_to(extend))
    | seq::to_std_array<3>();
}




//=============================================================================
face_godunov_data_t mhd_scheme_v2::godunov_fluxes(cell_primitive_variables_t pc, face_magnetic_flux_density_t bf)
{
    auto to_godunov_data = apply_to([] (auto ff, auto ef) -> mhd::godunov_data_t
    {
        return numeric::tuple_cat(ff, numeric::tuple(ef.component_1(), ef.component_2(), ef.component_3()));
    });
    auto [bf1, bf2, bf3] = bf;

    auto pf1 = nd::zip(pc | nd::select(0, 0, -1), pc | nd::select(0, 1), bf1);
    auto pf2 = nd::zip(pc | nd::select(1, 0, -1), pc | nd::select(1, 1), bf2);
    auto pf3 = nd::zip(pc | nd::select(2, 0, -1), pc | nd::select(2, 1), bf3);

    auto gf1 = pf1 | nd::map(rs1) | nd::map(to_godunov_data) | nd::to_shared();
    auto gf2 = pf2 | nd::map(rs2) | nd::map(to_godunov_data) | nd::to_shared();
    auto gf3 = pf3 | nd::map(rs3) | nd::map(to_godunov_data) | nd::to_shared();

    return std::array{gf1, gf2, gf3};
}




//=============================================================================
edge_electromotive_density_t mhd_scheme_v2::electromotive_forces(face_godunov_data_t gf)
{
    auto c1 = nd::map([] (auto g) { return numeric::get<5>(g); });
    auto c2 = nd::map([] (auto g) { return numeric::get<6>(g); });
    auto c3 = nd::map([] (auto g) { return numeric::get<7>(g); });
    auto a1 = nd::adjacent_mean(0);
    auto a2 = nd::adjacent_mean(1);
    auto a3 = nd::adjacent_mean(2);
    auto ev = nd::to_shared();
    auto [gf1, gf2, gf3] = gf;

    auto ee1a = gf2 | c1 | a3;
    auto ee2a = gf3 | c2 | a1;
    auto ee3a = gf1 | c3 | a2;
    auto ee1b = gf3 | c1 | a2;
    auto ee2b = gf1 | c2 | a3;
    auto ee3b = gf2 | c3 | a1;

    auto ee1 = 0.5 * (ee1a + ee1b);
    auto ee2 = 0.5 * (ee2a + ee2b);
    auto ee3 = 0.5 * (ee3a + ee3b);

    return std::array{ee1|ev, ee2|ev, ee3|ev};
}




//=============================================================================
cell_conserved_density_t mhd_scheme_v2::updated_conserved_density(
    cell_conserved_density_t uc,
    face_godunov_data_t gf,
    dimensional::unit_time dt,
    dimensional::unit_length dl)
{
    auto conserved_density_components = [] (auto g)
    {
        return numeric::tuple(
            numeric::get<0>(g),
            numeric::get<1>(g),
            numeric::get<2>(g),
            numeric::get<3>(g),
            numeric::get<4>(g));
    };

    auto [gf1, gf2, gf3] = gf;
    auto ff1 = gf1 | nd::map(conserved_density_components) | rt1;
    auto ff2 = gf2 | nd::map(conserved_density_components) | rt2;
    auto ff3 = gf3 | nd::map(conserved_density_components) | rt3;
    return (uc - mesh::divergence_difference(ff1, ff2, ff3) * dt / dl) | nd::to_shared();
}




//=============================================================================
face_magnetic_flux_density_t mhd_scheme_v2::curl(edge_electromotive_density_t ee, dimensional::unit_length dl)
{
    auto dt = dimensional::unit_time(1.0);
    auto [ee1, ee2, ee3] = ee;
    auto [cf1, cf2, cf3] = mesh::solenoidal_difference(ee1, ee2, ee3);

    return {
        cf1 * dt / dl | nd::to_shared(),
        cf2 * dt / dl | nd::to_shared(),
        cf3 * dt / dl | nd::to_shared(),
    };
}




//=============================================================================
face_magnetic_flux_density_t mhd_scheme_v2::updated_magnetic_flux_density(
    face_magnetic_flux_density_t bf,
    edge_electromotive_density_t ee,
    dimensional::unit_time dt,
    dimensional::unit_length dl)
{
    auto [bf1, bf2, bf3] = bf;
    auto [ee1, ee2, ee3] = ee;
    auto [cf1, cf2, cf3] = mesh::solenoidal_difference(ee1|rl1, ee2|rl2, ee3|rl3);

    return std::array{
        (bf1 - cf1 * dt / dl) | nd::to_shared(),
        (bf2 - cf2 * dt / dl) | nd::to_shared(),
        (bf3 - cf3 * dt / dl) | nd::to_shared()};
}
