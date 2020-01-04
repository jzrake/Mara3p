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




#include "mesh_sliding.hpp"
#include "scheme_plm_gradient.hpp"
#include "scheme_sedov2d.hpp"




//=============================================================================
dimensional::unit_area sedov::face_area(
    dimensional::unit_length r0, dimensional::unit_length r1,
    dimensional::unit_scalar t0, dimensional::unit_scalar t1)
{
    auto R0 = r0 * std::sin(t0);
    auto R1 = r1 * std::sin(t1);
    auto z0 = r0 * std::cos(t0);
    auto z1 = r1 * std::cos(t1);
    auto dR = R1 - R0;
    auto dz = z1 - z0;

    return M_PI * (R0 + R1) * dimensional::pow<1, 2>(dR * dR + dz * dz);
}

dimensional::unit_volume sedov::cell_volume(
    dimensional::unit_length r0, dimensional::unit_length r1,
    dimensional::unit_scalar t0, dimensional::unit_scalar t1)
{
    auto dcost = -(std::cos(t1) - std::cos(t0));
    return 2.0 * M_PI * (r1 * r1 * r1 - r0 * r0 * r0) / 3.0 * dcost;
}




//=============================================================================
sedov::radial_track_t sedov::generate_radial_track(
    dimensional::unit_length r0,
    dimensional::unit_length r1,
    dimensional::unit_scalar theta0, 
    dimensional::unit_scalar theta1)
{
    auto N = unsigned(M_PI / (theta1 - theta0));
    auto t = 0.5 * (theta0 + theta1);
    auto s = 0.2 * std::pow(std::cos(t * 10.0), 2.0);

    auto radii = nd::linspace(std::log10(r0.value), std::log10(r1.value), N)
    | nd::map([s, r0, r1] (double log10r)
    {
        auto r = dimensional::unit_length(std::pow(10.0, log10r));
        return r * (1.0 + (r / r0 - 1.0) * (r / r1 - 1.0) * s);
    });

    return {
        radii | nd::to_shared(),
        theta0,
        theta1,
    };
}




//=============================================================================
nd::shared_array<srhd::conserved_t, 1> sedov::generate_conserved(radial_track_t track)
{
    auto primitive = [] (dimensional::unit_length r, dimensional::unit_scalar t)
    {
        return srhd::primitive(1.0, 1.0);
    };
    auto rc = cell_center_radii(track);
    auto tc = cell_center_theta(track);
    auto dv = cell_volumes(track);

    return rc
    | nd::map(std::bind(primitive, std::placeholders::_1, tc))
    | nd::map(std::bind(srhd::conserved_density, std::placeholders::_1, 4. / 3))
    | nd::multiply(dv)
    | nd::to_shared();
}




//=============================================================================
nd::shared_array<srhd::primitive_t, 1> sedov::recover_primitive(
    radial_track_t track,
    nd::shared_array<srhd::conserved_t, 1> uc)
{
    auto p2c = std::bind(srhd::recover_primitive, std::placeholders::_1, 4. / 3, 0.0);
    auto dv = cell_volumes(track);
    return nd::to_shared((uc / dv) | nd::map(p2c));
}




//=============================================================================
nd::shared_array<sedov::primitive_per_length_t, 1> sedov::radial_gradient(
    radial_track_t track,
    nd::shared_array<srhd::primitive_t, 1> pc)
{
    auto plm = mara::plm_gradient(1.5);
    auto xc3 = cell_center_radii(track) | nd::adjacent_zip3();
    auto pc3 = pc | nd::adjacent_zip3();

    return nd::zip(xc3, pc3)
    | nd::map(plm)
    | nd::extend_zeros()
    | nd::to_shared();
}




//=============================================================================
nd::shared_array<sedov::radial_godunov_data_t, 1> sedov::radial_godunov_data(
        radial_track_t track,
        nd::shared_array<srhd::primitive_t, 1> pc,
        nd::shared_array<primitive_per_length_t, 1> dc)
{
    auto nhat = geometric::unit_vector_on(1);
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_across_contact_t();
    // auto mode = srhd::riemann_solver_mode_hllc_fluxes_t();

    auto riemann = [nhat, mode] (srhd::primitive_t pl, srhd::primitive_t pr) -> radial_godunov_data_t
    {
        return srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode);
        // return std::pair(srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode), dimensional::unit_velocity(0.0));
    };

    return pc
    | nd::adjacent_zip()
    | nd::map(util::apply_to(riemann))
    | nd::extend_zero_gradient()
    | nd::to_shared();
}




//=============================================================================
nd::shared_array<sedov::polar_godunov_data_t, 1> sedov::polar_godunov_data(track_data_t L, track_data_t R)
{
    auto tl = std::get<0>(L);
    auto tr = std::get<0>(R);
    auto pl = std::get<1>(L);
    auto pr = std::get<1>(R);

    auto nhat = geometric::unit_vector_on(2);
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_t();
    auto face = polar_faces(tl, tr);

    return nd::make_array(nd::indexing([tl, tr, pl, pr, face, nhat, mode] (nd::uint i) -> polar_godunov_data_t
    {
        if (face(i).il && face(i).ir)
        {
            auto il = face(i).il.value();
            auto ir = face(i).ir.value();
            auto ff = srhd::riemann_solver(pl(il), pr(ir), nhat, 4. / 3, mode);
            auto da = face_area(face(i).trailing, face(i).leading, tl.theta1, tr.theta0);
            return {ff, da, il, ir};
        }
        return {srhd::flux_vector_t(), dimensional::unit_area(0.0), 0, 0};
    }), shape(face))
    | nd::to_shared();
}




//=============================================================================
nd::shared_array<srhd::conserved_t, 1> sedov::delta_conserved(radial_track_t track,
    nd::shared_array<srhd::primitive_t, 1> pc,
    nd::shared_array<radial_godunov_data_t, 1> ff,
    nd::shared_array<polar_godunov_data_t, 1> gfl,
    nd::shared_array<polar_godunov_data_t, 1> gfr,
    dimensional::unit_time dt)
{
    auto [fhat, vhat] = nd::unzip(ff);
    auto da = radial_face_areas(track);

    auto df = (fhat * da) | nd::adjacent_diff();
    auto [ghat_l, da_l, il_l, ir_l] = nd::unzip(gfl);
    auto [ghat_r, da_r, il_r, ir_r] = nd::unzip(gfr);
    auto gl = mesh::bin_values(ghat_l * da_l, ir_l, size(df));
    auto gr = mesh::bin_values(ghat_r * da_r, il_r, size(df));
    auto dg = gr - gl;

    auto source_terms = util::apply_to(std::bind(
        srhd::spherical_geometry_source_terms,
        std::placeholders::_1,
        std::placeholders::_2,
        cell_center_theta(track),
        4. / 3));

    auto sc = nd::zip(pc, cell_center_radii(track))
    | nd::map(source_terms)
    | nd::multiply(cell_volumes(track));

    return nd::to_shared((-(df + dg) + sc) * dt);
}




//=============================================================================
nd::shared_array<dimensional::unit_length, 1> sedov::delta_face_positions(
    nd::shared_array<radial_godunov_data_t, 1> ff,
    dimensional::unit_time dt)
{
    auto [fhat, vhat] = nd::unzip(ff);
    return nd::to_shared(vhat * dt);
}




//=============================================================================
std::pair<sedov::radial_track_t, nd::shared_array<srhd::primitive_t, 1>> sedov::extend(
    radial_track_t track,
    nd::shared_array<srhd::primitive_t, 1> pc)
{
    return {
        {
            track.face_radii | nd::extend_extrap() | nd::to_shared(),
            track.theta0,
            track.theta1,
        },
        pc | nd::extend_uniform(srhd::primitive(1.0, 1.0)) | nd::to_shared(),
    };
}
