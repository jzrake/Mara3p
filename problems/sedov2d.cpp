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




#include <iostream>
#include "app_vtk.hpp"
#include "core_dimensional.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"
#include "mesh_sliding.hpp"
#include "physics_srhd.hpp"
#include "scheme_plm_gradient.hpp"




//=============================================================================
namespace sedov {




//=============================================================================
struct solution_state_t
{
    rational::number_t                            iteration;
    dimensional::unit_time                        time;
    nd::shared_array<dimensional::unit_scalar, 1> polar_vertices;
};

struct radial_track_t
{
    nd::shared_array<dimensional::unit_length, 1> face_radii;
    dimensional::unit_scalar                      theta0;
    dimensional::unit_scalar                      theta1;
};

using primitive_per_length_t = decltype(srhd::primitive_t() / dimensional::unit_length());
using primitive_per_scalar_t = decltype(srhd::primitive_t() / dimensional::unit_scalar());
using radial_godunov_data_t  = std::tuple<srhd::flux_vector_t, dimensional::unit_velocity>;
using polar_godunov_data_t   = std::tuple<srhd::flux_vector_t, dimensional::unit_area, nd::uint, nd::uint>;
using track_data_t           = std::tuple<radial_track_t, nd::shared_array<srhd::primitive_t, 1>, nd::shared_array<primitive_per_length_t, 1>>;




//=============================================================================
dimensional::unit_area face_area(
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

dimensional::unit_volume cell_volume(
    dimensional::unit_length r0, dimensional::unit_length r1,
    dimensional::unit_scalar t0, dimensional::unit_scalar t1)
{
    auto dcost = -(std::cos(t1) - std::cos(t0));
    return 2.0 * M_PI * (r1 * r1 * r1 - r0 * r0 * r0) / 3.0 * dcost;
}




//=============================================================================
auto radial_face_areas(radial_track_t track)
{
    return track.face_radii
    | nd::map([t0=track.theta0, t1=track.theta1] (auto r)
    {
        return face_area(r, r, t0, t1);
    });
}

auto cell_volumes(radial_track_t track)
{
    return track.face_radii
    | nd::adjacent_zip()
    | nd::map(util::apply_to([t0=track.theta0, t1=track.theta1] (auto r0, auto r1)
    {
        return cell_volume(r0, r1, t0, t1);
    }));
}

auto cell_center_radii(radial_track_t track)
{
    return track.face_radii | nd::adjacent_mean();
}

auto cell_center_theta(radial_track_t track)
{
    return 0.5 * (track.theta0 + track.theta1);
}

auto polar_faces(radial_track_t L, radial_track_t R)
{
    return mesh::transverse_faces(L.face_radii, R.face_radii);
}




//=============================================================================
nd::shared_array<radial_track_t, 1> initial_radial_tracks(
    unsigned track_count,
    dimensional::unit_length r0,
    dimensional::unit_length r1)
{
    auto theta = nd::linspace(0.0, M_PI, track_count + 1);

    return nd::range(track_count) | nd::map([r0, r1, track_count, theta] (nd::uint j) -> radial_track_t
    {
        auto N = unsigned(std::log10(r1 / r0)) * track_count;
        auto s = j % 2 ? 1.0 : 1.0 + 1.0 / track_count;

        auto radii = nd::linspace(std::log10(s * r0.value), std::log10(s * r1.value), N)
        | nd::map([] (double log10r) { return std::pow(10.0, log10r); })
        | nd::construct<dimensional::unit_length>();

        return {
            radii | nd::to_shared(),
            theta(j + 0),
            theta(j + 1),
        };
    }) | nd::to_shared();
}

nd::shared_array<srhd::primitive_t, 1> initial_primitive(radial_track_t track)
{
    auto primitive = [] (dimensional::unit_length r, dimensional::unit_scalar t)
    {
        return srhd::primitive(1.0, 1.0);
    };
    auto rc = cell_center_radii(track);
    auto tc = cell_center_theta(track);

    return rc | nd::map(std::bind(primitive, std::placeholders::_1, tc)) | nd::to_shared();
}




//=============================================================================
nd::shared_array<srhd::conserved_t, 1> conserved(radial_track_t track, nd::shared_array<srhd::primitive_t, 1> pc)
{
    auto c2p = std::bind(srhd::conserved_density, std::placeholders::_1, 4. / 3);
    auto dv = cell_volumes(track);
    return (map(pc, c2p) * dv) | nd::to_shared();
}

nd::shared_array<srhd::primitive_t, 1> recover_primitive(radial_track_t track, nd::shared_array<srhd::conserved_t, 1> uc)
{
    auto p2c = std::bind(srhd::recover_primitive, std::placeholders::_1, 4. / 3, 0.0);
    auto dv = cell_volumes(track);
    return (uc / dv) | nd::map(p2c) | nd::to_shared();
}

nd::shared_array<primitive_per_length_t, 1> radial_gradient(radial_track_t track, nd::shared_array<srhd::primitive_t, 1> pc)
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
nd::shared_array<radial_godunov_data_t, 1> radial_godunov_data(
        radial_track_t track,
        nd::shared_array<srhd::primitive_t, 1> pc,
        nd::shared_array<primitive_per_length_t, 1> dc)
{
    auto nhat = geometric::unit_vector_on(1);
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_across_contact_t();

    auto riemann = [nhat, mode] (srhd::primitive_t pl, srhd::primitive_t pr) -> radial_godunov_data_t
    {
        return srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode);
    };

    return pc
    | nd::adjacent_zip()
    | nd::map(util::apply_to(riemann))
    | nd::to_shared();
}




//=============================================================================
nd::shared_array<polar_godunov_data_t, 1> polar_godunov_data(track_data_t L, track_data_t R)
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
        return {{}, {}, 0, 0};
    }), shape(face))
    | nd::to_shared();
}




//=============================================================================
auto quad_mesh(nd::shared_array<radial_track_t, 1> tracks)
{
    using vec3d = geometric::euclidean_vector_t<dimensional::unit_length>;
    auto total_cells = tracks | nd::map([] (auto track) { return size(track.face_radii) - 1; }) | nd::sum();

    auto vertices = nd::make_unique_array<vec3d>(nd::uivec(total_cells * 4));
    auto n = nd::uint(0);

    for (std::size_t j = 0; j < size(tracks); ++j)
    {
        for (std::size_t i = 0; i < size(tracks(j).face_radii) - 1; ++i)
        {
            auto t0 = tracks(j).theta0;
            auto t1 = tracks(j).theta1;
            auto r0 = tracks(j).face_radii(i + 0);
            auto r1 = tracks(j).face_radii(i + 1);

            vertices(n++) = vec3d{r0 * std::sin(t0), 0.0, r0 * std::cos(t0)};
            vertices(n++) = vec3d{r0 * std::sin(t1), 0.0, r0 * std::cos(t1)};
            vertices(n++) = vec3d{r1 * std::sin(t1), 0.0, r1 * std::cos(t1)};
            vertices(n++) = vec3d{r1 * std::sin(t0), 0.0, r1 * std::cos(t0)};
        }
    }

    auto indexes = nd::range(total_cells)
    | nd::map([] (int n)
    {
        return std::array{4 * n + 0, 4 * n + 1, 4 * n + 2, 4 * n + 3};
    });
    return std::pair(nd::make_shared_array(std::move(vertices)), nd::to_shared(indexes));
}

} // namespace sedov




//=============================================================================
int main()
{
    auto num_tracks = 100;

    auto tracks = sedov::initial_radial_tracks(num_tracks, 1.0, 10.0);
    auto face_areas   = tracks | nd::map(sedov::radial_face_areas);
    auto cell_volumes = tracks | nd::map(sedov::cell_volumes);

    auto [vertices, indexes] = sedov::quad_mesh(tracks);

    auto rc = tracks | nd::map(sedov::cell_center_radii) | nd::flat();
    auto pc = tracks | nd::map(sedov::initial_primitive);

    auto uc = nd::zip(tracks, pc) | nd::map(util::apply_to(sedov::conserved));
    auto dc = nd::zip(tracks, pc) | nd::map(util::apply_to(sedov::radial_gradient));
    auto gf = nd::zip(tracks, pc, dc) | nd::map(util::apply_to(sedov::radial_godunov_data));
    auto hf = nd::zip(tracks, pc, dc) | nd::adjacent_zip() | nd::map(util::apply_to(sedov::polar_godunov_data));

    auto [fhat, da, il, ir] = nd::unzip(hf(0));
    auto l0 = mesh::bin_values(-fhat * da, il, size(uc(0)));
    auto l1 = mesh::bin_values(+fhat * da, ir, size(uc(1)));

    // vtk::write(std::cout, "Grid", vertices, indexes, std::pair("cell radius", rc));

    return 0;
}
