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
#include <fstream>
#include "app_config.hpp"
#include "app_vtk.hpp"
#include "core_util.hpp"
#include "scheme_sedov2d.hpp"




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("nr",                   256)   // number of radial zones, per decade
    .item("tfinal",              10.0)   // time to stop the simulation
    .item("print",                 10)   // the number of iterations between terminal outputs
    .item("dfi",                  1.5)   // output interval (constant multiplier)
    .item("rk_order",               2)   // Runge-Kutta order (1, 2, or 3)
    .item("cfl",                  0.5)   // courant number
    .item("plm_theta",            1.5)   // PLM parameter
    .item("router",               1e1)   // outer boundary radius
    .item("move",                   1)   // whether to move the cells
    .item("envelop_mdot_index",   4.0)   // alpha in envelop mdot(t) = (t / t0)^alpha
    .item("envelop_u_index",     0.22)   // psi in envelop u(m) = u1 (m / m1)^(-psi)
    .item("envelop_u",          10.00)   // maximum gamma-beta in outer envelop
    .item("engine_mdot",          1e4)   // engine mass rate
    .item("engine_onset",        50.0)   // the engine onset time [inner boundary light-crossing time]
    .item("engine_duration",    100.0)   // the engine duration   [inner boundary light-crossing time]
    .item("engine_u",            10.0);  // the engine gamma-beta
}




//=============================================================================
struct solution_state_t
{
    rational::number_t                                          iteration;
    dimensional::unit_time                                      time;
    nd::shared_array<sedov::radial_track_t, 1>                  tracks;
    nd::shared_array<nd::shared_array<srhd::conserved_t, 1>, 1> conserved;
};




srhd::primitive_t initial_primitive(dimensional::unit_length r, dimensional::unit_scalar q, dimensional::unit_time t)
{
    auto r0 = dimensional::unit_length(1.0);
    auto d0 = dimensional::unit_mass_density(std::pow(r / r0, -2.0));
    auto c2 = srhd::light_speed * srhd::light_speed;
    auto ur = 0.1 + 0.2 * (std::pow(std::cos(q), 2) + std::pow(std::cos(M_PI - q), 2));
    return srhd::primitive(d0, ur, 0.0, 0.0, 1e-6 * d0 * c2);
}




//=============================================================================
solution_state_t initial_solution(const mara::config_t& cfg)
{
    using namespace std::placeholders;

    auto r0 = dimensional::unit_length(1.0);
    auto r1 = r0 * cfg.get_double("router");
    auto p0 = [] (dimensional::unit_length r, dimensional::unit_scalar q) -> srhd::primitive_t
    {
        return initial_primitive(r, q, 0.0);
    };

    auto num_tracks = cfg.get_int("nr");
    auto tr = nd::linspace(0.0, M_PI, num_tracks + 1)
    | nd::adjacent_zip()
    | nd::map(util::apply_to(std::bind(sedov::generate_radial_track, r0, r1, _1, _2)));
    auto u0 = tr | nd::map(std::bind(sedov::generate_conserved, _1, p0));

    return {
        0,
        0.0,
        nd::to_shared(tr),
        nd::to_shared(u0),
    };
}

sedov::radial_godunov_data_t inner_bc(sedov::radial_track_t track, nd::shared_array<srhd::primitive_t, 1> pc, dimensional::unit_time t)
{
    auto vel = dimensional::unit_velocity(0.0);
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_moving_face_t{vel};
    auto nhat = geometric::unit_vector_on(1);
    auto pl = initial_primitive(front(track.face_radii), cell_center_theta(track), t);
    auto pr = front(pc);
    return std::pair(srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode), vel);
}

sedov::radial_godunov_data_t outer_bc(sedov::radial_track_t track, nd::shared_array<srhd::primitive_t, 1> pc, dimensional::unit_time t)
{
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_across_contact_t();
    auto nhat = geometric::unit_vector_on(1);
    auto pl = back(pc);
    auto pr = initial_primitive(back(track.face_radii), cell_center_theta(track), t);
    return srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode);
}

solution_state_t remesh(solution_state_t solution)
{
    using namespace std::placeholders;

    auto t0 = solution.tracks;
    auto uc = solution.conserved;

    auto [t1, u1] = nd::unzip(nd::zip(t0, uc) | nd::map(util::apply_to(std::bind(sedov::refine, _1, _2, dimensional::unit_scalar(1.33)))));

    return {solution.iteration, solution.time, t1|nd::to_shared(), u1|nd::to_shared()};
}




//=============================================================================
solution_state_t advance(solution_state_t solution, dimensional::unit_time dt)
{
    using namespace std::placeholders;

    auto t0 = solution.tracks;
    auto uc = solution.conserved;
    auto pc = nd::zip(t0, uc)             | nd::map(util::apply_to(sedov::recover_primitive));
    auto fi = nd::zip(t0, pc)             | nd::map(util::apply_to(std::bind(inner_bc, _1, _2, solution.time)));
    auto fo = nd::zip(t0, pc)             | nd::map(util::apply_to(std::bind(outer_bc, _1, _2, solution.time)));
    auto dc = nd::zip(t0, pc)             | nd::map(util::apply_to(sedov::radial_gradient));
    auto ff = nd::zip(t0, pc, dc, fi, fo) | nd::map(util::apply_to(sedov::radial_godunov_data));
    auto dr = ff                          | nd::map(std::bind(sedov::delta_face_positions, _1, dt));

    auto gf = nd::zip(t0, pc, dc)
    | nd::adjacent_zip()
    | nd::map(util::apply_to(sedov::polar_godunov_data))
    | nd::extend_uniform(nd::shared_array<sedov::polar_godunov_data_t, 1>{});

    auto du = nd::range(size(uc))
    | nd::map([t0, pc, ff, gf, dt] (nd::uint j)
    {
        return sedov::delta_conserved(t0(j), pc(j), ff(j), gf(j), gf(j + 1), dt);
    });

    auto t1 = nd::zip(t0, dr) | nd::map(util::apply_to([] (auto track, auto dr)
    {
        return sedov::radial_track_t{
            nd::to_shared(track.face_radii + dr),
            track.theta0,
            track.theta1,
        };
    }));

    auto u1 = nd::zip(uc, du) | nd::map(util::apply_to([] (auto uc, auto du)
    {
        return nd::to_shared(uc + du);
    }));

    return remesh({
        solution.iteration + 1,
        solution.time + dt,
        t1 | nd::to_shared(),
        u1 | nd::to_shared(),
    });
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
int main(int argc, const char* argv[])
{
    auto cfg = config_template()
    .create()
    .update(mara::argv_to_string_map(argc, argv));

    mara::pretty_print(std::cout, "config", cfg);

    auto solution = initial_solution(cfg);

    auto tfinal = dimensional::unit_time  (cfg.get_double("tfinal"));
    auto cfl    = dimensional::unit_scalar(cfg.get_double("cfl"));
    auto dt     = cfl * nd::min(solution.tracks | nd::map(sedov::minimum_spacing)) / srhd::light_speed;

    auto vtk_count = 0;
    output_vtk(solution, vtk_count);

    while (solution.time < tfinal)
    {
        std::printf("[%06lu] t=%.4f\n", long(solution.iteration), solution.time.value);
        solution = advance(solution, dt);

        if (long(solution.iteration) % 1 == 0)
        {
            output_vtk(solution, vtk_count++);            
        }
    }
    return 0;
}
