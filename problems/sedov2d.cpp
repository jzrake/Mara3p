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




#include <fstream>
#include "app_vtk.hpp"
#include "core_util.hpp"
#include "scheme_sedov2d.hpp"




//=============================================================================
struct solution_state_t
{
    rational::number_t                                          iteration;
    dimensional::unit_time                                      time;
    nd::shared_array<sedov::radial_track_t, 1>                  tracks;
    nd::shared_array<nd::shared_array<srhd::conserved_t, 1>, 1> conserved;
};




//=============================================================================
solution_state_t initial_solution()
{
    auto num_tracks = 20;
    auto t0 = sedov::generate_radial_tracks(num_tracks, 1.0, 10.0);
    auto u0 = t0 | nd::map(sedov::generate_conserved) | nd::to_shared();

    return {
        0,
        0.0,
        t0,
        u0,
    };
}




//=============================================================================
solution_state_t advance(solution_state_t solution, dimensional::unit_time dt)
{
    auto t0 = solution.tracks;
    auto u0 = solution.conserved;

    auto p0 = nd::zip(t0, u0) | nd::map(util::apply_to(sedov::recover_primitive));
    auto [te, pc] = nd::unzip(nd::zip(t0, p0) | nd::map(util::apply_to(sedov::extend)));
    auto dc = nd::zip(te, pc) | nd::map(util::apply_to(sedov::radial_gradient));
    auto ff = nd::zip(te, pc, dc) | nd::map(util::apply_to(sedov::radial_godunov_data));
    auto gf = nd::zip(te, pc, dc)
    | nd::adjacent_zip()
    | nd::map(util::apply_to(sedov::polar_godunov_data))
    | nd::extend_uniform(nd::shared_array<sedov::polar_godunov_data_t, 1>{})
    | nd::to_shared();

    auto dr = ff | nd::map(std::bind(sedov::delta_face_positions, std::placeholders::_1, dt));
    auto du = nd::range(size(u0))
    | nd::map([t0, p0, ff, gf, dt] (nd::uint j)
    {
        return sedov::delta_conserved(t0(j), p0(j), ff(j), gf(j), gf(j + 1), dt);
    });

    auto t1 = nd::zip(t0, dr) | nd::map(util::apply_to([] (auto track, auto dr)
    {
        return sedov::radial_track_t{
            nd::to_shared(track.face_radii + dr),
            track.theta0,
            track.theta1,
        };
    })) | nd::to_shared();

    auto u1 = nd::zip(u0, du) | nd::map(util::apply_to([] (auto u0, auto du)
    {
        return nd::to_shared(u0 + du);
    })) | nd::to_shared();

    return {
        solution.iteration + 1,
        solution.time + dt,
        t1,
        u1,
    };
}




//=============================================================================
auto quad_mesh(nd::shared_array<sedov::radial_track_t, 1> tracks)
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




//=============================================================================
void output_vtk(solution_state_t solution, unsigned count)
{
    auto fname = util::format("primitive.%04u.vtk", count);
    auto outf = std::ofstream(fname, std::ios_base::out);
    auto [vertices, indexes] = quad_mesh(solution.tracks);

    auto p0 = nd::zip(solution.tracks, solution.conserved) | nd::map(util::apply_to(sedov::recover_primitive));
    auto density = p0 | nd::flat() | nd::map(srhd::mass_density) | nd::to_shared();
    auto ur      = p0 | nd::flat() | nd::map(srhd::gamma_beta_1) | nd::to_shared();

    std::printf("output %s\n", fname.data());

    vtk::write(outf, "Grid", vertices, indexes,
        std::pair("density", density),
        std::pair("radial-gamma-beta", ur));
}




//=============================================================================
int main()
{
    auto dt = dimensional::unit_time(0.001);
    auto solution = initial_solution();

    output_vtk(solution, 0);

    while (long(solution.iteration) < 10)
    {
        std::printf("[%06lu] t=%.4f\n", long(solution.iteration), solution.time.value);
        solution = advance(solution, dt);
    }
    output_vtk(solution, 1);

    return 0;
}
