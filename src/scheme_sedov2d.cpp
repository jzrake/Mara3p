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

static const double temperature_floor = 1e-6;




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
dimensional::unit_length sedov::minimum_spacing(radial_track_t track)
{
    return nd::min(track.face_radii | nd::adjacent_diff());
}

dimensional::unit_scalar sedov::cell_center_theta(radial_track_t track)
{
    return 0.5 * (track.theta0 + track.theta1);
}

nd::shared_array<dimensional::unit_length, 1> sedov::cell_center_radii(radial_track_t track)
{
    return track.face_radii
    | nd::adjacent_mean()
    | nd::to_shared();
}

nd::shared_array<dimensional::unit_area, 1> sedov::radial_face_areas(radial_track_t track)
{
    return track.face_radii
    | nd::map([t0=track.theta0, t1=track.theta1] (auto r)
    {
        return face_area(r, r, t0, t1);
    })
    | nd::to_shared();
}

nd::shared_array<dimensional::unit_volume, 1> sedov::cell_volumes(radial_track_t track)
{
    return track.face_radii
    | nd::adjacent_zip()
    | nd::map(util::apply_to([t0=track.theta0, t1=track.theta1] (auto r0, auto r1)
    {
        return cell_volume(r0, r1, t0, t1);
    }))
    | nd::to_shared();
}




//=============================================================================
sedov::radial_track_t sedov::generate_radial_track(
    dimensional::unit_length r0,
    dimensional::unit_length r1,
    dimensional::unit_scalar theta0,
    dimensional::unit_scalar theta1,
    dimensional::unit_scalar aspect)
{
    auto N = unsigned(M_PI / (theta1 - theta0) / aspect);

    auto radii = nd::linspace(std::log10(r0.value), std::log10(r1.value), N)
    | nd::map([] (double log10r)
    {
        return dimensional::unit_length(std::pow(10.0, log10r));
    });

    return {
        radii | nd::to_shared(),
        theta0,
        theta1,
    };
}




//=============================================================================
nd::shared_array<srhd::conserved_t, 1> sedov::generate_conserved(radial_track_t track, primitive_function_t primitive)
{
    auto rc = cell_center_radii(track);
    auto qc = cell_center_theta(track);
    auto dv = cell_volumes(track);

    return rc
    | nd::map(std::bind(primitive, std::placeholders::_1, qc))
    | nd::map(std::bind(srhd::conserved_density, std::placeholders::_1, 4. / 3))
    | nd::multiply(dv)
    | nd::to_shared();
}




//=============================================================================
nd::shared_array<srhd::primitive_t, 1> sedov::recover_primitive(
    radial_track_t track,
    nd::shared_array<srhd::conserved_t, 1> uc)
{
    auto p2c = std::bind(srhd::recover_primitive, std::placeholders::_1, 4. / 3, temperature_floor);
    auto dv = cell_volumes(track);

    try {
        return nd::to_shared((uc / dv) | nd::map(p2c));
    }
    catch (const std::exception& e) {
        nd::uint index = 0;

        for (auto [i, qc] : nd::zip(nd::range(size(uc)), uc / dv))
        {
            try {
                p2c(qc);
            }
            catch (...) {
                index = i;
                break;
            }
        }
        throw std::runtime_error(util::format("%s at index %lu / %lu r = [%lf, %lf] theta = [%lf, %lf] ",
            e.what(),
            index,
            size(uc),
            track.face_radii(index).value,
            track.face_radii(index + 1).value,
            track.theta0.value,
            track.theta1.value));
    }
}




//=============================================================================
nd::shared_array<sedov::primitive_per_length_t, 1> sedov::radial_gradient(
    radial_track_t track,
    nd::shared_array<srhd::primitive_t, 1> pc,
    bool use_plm)
{
    auto plm = mara::plm_gradient(use_plm ? 1.0 : 0.0);
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
        nd::shared_array<primitive_per_length_t, 1> dc,
        radial_godunov_data_t inner_boundary_data,
        radial_godunov_data_t outer_boundary_data)
{
    auto nhat = geometric::unit_vector_on(1);
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_across_contact_t();

    auto riemann = [nhat, mode] (srhd::primitive_t pl, srhd::primitive_t pr) -> radial_godunov_data_t
    {
        return srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode);
    };

    auto dr = track.face_radii | nd::adjacent_diff();
    auto pl = select(pc + 0.5 * dc * dr, 0, 0, -1);
    auto pr = select(pc - 0.5 * dc * dr, 0, 1);

    return nd::zip(pl, pr)
    | nd::map(util::apply_to(riemann))
    | nd::extend_uniform_lower(inner_boundary_data)
    | nd::extend_uniform_upper(outer_boundary_data)
    | nd::to_shared();
}




//=============================================================================
srhd::primitive_t sedov::sample(track_data_t track_data, dimensional::unit_length r, nd::uint index, srhd::primitive_t fallback)
{
    auto [tr, pc, dc] = track_data;

    if (r <= front(tr.face_radii) || r >= back(tr.face_radii))
    {
        return fallback;
    }
    index = std::min(index, size(pc));

    while (r < tr.face_radii(index + 0)) --index;
    while (r > tr.face_radii(index + 1)) ++index;

    auto r0 = tr.face_radii(index + 0);
    auto r1 = tr.face_radii(index + 1);
    auto rc = 0.5 * (r0 + r1);

    return pc(index) + dc(index) * (r - rc);
}




//=============================================================================
nd::shared_array<sedov::polar_godunov_data_t, 1> sedov::polar_godunov_data(track_data_t t0, track_data_t t1, track_data_t t2, track_data_t t3, bool use_plm)
{
    auto nhat = geometric::unit_vector_on(2);
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_t();
    auto face = mesh::transverse_faces(std::get<0>(t1).face_radii, std::get<0>(t2).face_radii);
    auto plm  = mara::plm_gradient(1.0);

    return nd::make_array(nd::indexing([=] (nd::uint i) -> polar_godunov_data_t
    {
        if (face(i).il && face(i).ir)
        {
            auto il = face(i).il.value();
            auto ir = face(i).ir.value();
            auto ri = face(i).trailing;
            auto ro = face(i).leading;
            auto ql = std::get<0>(t1).theta1;
            auto qr = std::get<0>(t2).theta0; // NOTE: ql and qr must be equal
            auto rf = 0.5 * (ri + ro);

            if (! use_plm || size(std::get<0>(t0).face_radii) == 0 || size(std::get<0>(t3).face_radii) == 0)
            {
                // If either the left-most or right-most track data is missing, then
                // forego extrapolation in the polar direction.
                auto pl = sample(t1, rf, il, {});
                auto pr = sample(t2, rf, ir, {});
                auto ff = srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode);
                auto da = face_area(ri, ro, ql, qr);

                return {ff, da, il, ir};
            }

            auto q0 = cell_center_theta(std::get<0>(t0));
            auto q1 = cell_center_theta(std::get<0>(t1));
            auto q2 = cell_center_theta(std::get<0>(t2));
            auto q3 = cell_center_theta(std::get<0>(t3));
            auto p1 = sample(t1, rf, il, {});
            auto p2 = sample(t2, rf, ir, {});
            auto p0 = sample(t0, rf, il, p1);
            auto p3 = sample(t3, rf, ir, p2);
            auto cl = plm(std::tuple(std::tuple(q0, q1, q2), std::tuple(p0, p1, p2)));
            auto cr = plm(std::tuple(std::tuple(q1, q2, q3), std::tuple(p1, p2, p3)));
            auto pl = p1 + (ql - q1) * cl;
            auto pr = p2 + (qr - q2) * cr;
            auto ff = srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode);
            auto da = face_area(ri, ro, ql, qr);

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
std::pair<sedov::radial_track_t, nd::shared_array<srhd::conserved_t, 1>> sedov::remesh(
    radial_track_t track,
    nd::shared_array<srhd::conserved_t, 1> uc,
    dimensional::unit_scalar maximum_cell_aspect_ratio,
    dimensional::unit_scalar minimum_cell_aspect_ratio,
    dimensional::unit_length inner_boundary_radius)
{
    auto rf = track.face_radii;
    auto rc = rf | nd::adjacent_mean();
    auto dr = rf | nd::adjacent_diff();
    auto ds = rc * (track.theta1 - track.theta0);
    auto aspect = dr / ds;
    auto imin = nd::argmin(aspect)[0];
    auto imax = nd::argmax(aspect)[0];

    auto make_track = [track, minimum_cell_aspect_ratio, maximum_cell_aspect_ratio, inner_boundary_radius] (auto rf, auto uc)
    {
        return remesh({rf, track.theta0, track.theta1}, uc, minimum_cell_aspect_ratio, maximum_cell_aspect_ratio, inner_boundary_radius);
    };

    auto make_track_and_return = [track] (auto rf, auto uc)
    {
        return std::make_pair(radial_track_t{rf, track.theta0, track.theta1}, uc);
    };

    if (front(rf) < inner_boundary_radius)
    {
        return make_track(rf | nd::select(0, 1) | nd::to_shared(), uc | nd::select(0, 1) | nd::to_shared());
    }
    if (aspect(imax) > maximum_cell_aspect_ratio)
    {
        return std::apply(make_track_and_return, nd::add_partition(rf, uc, imax));
    }
    if (aspect(imin) < minimum_cell_aspect_ratio)
    {
        if (imin == 0)
        {
            return std::apply(make_track, nd::remove_partition(rf, uc, 1));
        }
        if (imin + 1 == size(aspect))
        {
            return std::apply(make_track, nd::remove_partition(rf, uc, imin));
        }
        if (aspect(imin - 1) <= aspect(imin + 1))
        {
            return std::apply(make_track, nd::remove_partition(rf, uc, imin));
        }
        if (aspect(imin - 1) >= aspect(imin + 1))
        {
            return std::apply(make_track, nd::remove_partition(rf, uc, imin + 1));
        }
    }
    return std::pair(track, uc);
}




//=============================================================================
sedov::track_data_t sedov::copy_track_data(sedov::track_data_t tr)
{
    return sedov::track_data_t{
        {
            std::get<0>(tr).face_radii | nd::to_shared(),
            std::get<0>(tr).theta0,
            std::get<0>(tr).theta1
        },
        std::get<1>(tr) | nd::to_shared(),
        std::get<2>(tr) | nd::to_shared(),
    };
}




//=============================================================================
std::tuple<sedov::track_data_t, sedov::track_data_t, sedov::track_data_t, sedov::track_data_t>
sedov::copy_track_data4(std::tuple<sedov::track_data_t, sedov::track_data_t, sedov::track_data_t, sedov::track_data_t> trs)
{
    return std::tuple(
        copy_track_data(std::get<0>(trs)),
        copy_track_data(std::get<1>(trs)),
        copy_track_data(std::get<2>(trs)),
        copy_track_data(std::get<3>(trs))
    );
}
