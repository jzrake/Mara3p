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
#include "app_control.hpp"
#include "app_vtk.hpp"
#include "core_util.hpp"
#include "model_wind.hpp"
#include "scheme_sedov2d.hpp"

using namespace dimensional;
using namespace std::placeholders;




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("nt",                   256)   // number of radial tracks
    .item("rk",                   2)     // Runge-Kutta order (1, 2, or 3)
    .item("aspect",               4.0)   // aspect ratio of the longest cells allowed
    .item("focus",                0.0)   // amount by which to increase resolution near the poles [0, 1)
    .item("tfinal",              10.0)   // time to stop the simulation
    // .item("print",                 10)   // the number of iterations between terminal outputs
    .item("vtk",                   50)   // number of iterations between VTK outputs
    // .item("dfi",                  1.5)   // output interval (constant multiplier)
    // .item("rk_order",               2)   // Runge-Kutta order (1, 2, or 3)
    .item("cfl",                  0.5)   // courant number
    // .item("plm_theta",            1.5)   // PLM parameter
    .item("router",               1e1)   // outer boundary radius
    .item("envelop_mdot_index",   4.0)   // alpha in envelop mdot(t) = (t / t0)^alpha
    .item("envelop_u_index",     0.22)   // psi in envelop u(m) = u1 (m / m1)^(-psi)
    .item("envelop_u",          10.00)   // maximum gamma-beta in outer envelop
    .item("engine_mdot",          1e4)   // engine mass rate
    .item("engine_onset",        50.0)   // the engine onset time [inner boundary light-crossing time]
    .item("engine_duration",    100.0)   // the engine duration   [inner boundary light-crossing time]
    .item("engine_u",            10.0);  // the engine gamma-beta
}




//=============================================================================
struct solution_t
{
    rational::number_t                                          iteration;
    dimensional::unit_time                                      time;
    nd::shared_array<sedov::radial_track_t, 1>                  tracks;
    nd::shared_array<nd::shared_array<srhd::conserved_t, 1>, 1> conserved;
};

solution_t weighted_sum(solution_t s, solution_t t, rational::number_t b)
{
    auto result = solution_t();

    result.iteration = s.iteration  *        b  + t.iteration *       (1 - b);
    result.time      = s.time       * double(b) + t.time      * double(1 - b);
    result.tracks    = nd::zip(s.tracks, t.tracks)
    | nd::map(util::apply_to([b] (auto s, auto t) -> sedov::radial_track_t
    {
        return {
            nd::to_shared(s.face_radii * double(b) + t.face_radii * double(1 - b)),
            s.theta0,
            s.theta1
        };
    })) | nd::to_shared();
    result.conserved = nd::zip(s.conserved, t.conserved)
    | nd::map(util::apply_to([b] (auto s, auto t)
    {
        return nd::to_shared(s * double(b) + t * double(1 - b));
    })) | nd::to_shared();

    return result;
}




//=============================================================================
auto wind_mass_loss_rate(const mara::config_t& cfg, unit_scalar)
{
    auto envelop_mdot = unit_mass_rate(1.0);
    auto engine_mdot  = unit_mass_rate(cfg.get_double("engine_mdot"));
    auto engine_onset = unit_time(cfg.get_double("engine_onset"));
    auto alpha        = cfg.get_double("envelop_mdot_index");

    return [envelop_mdot, engine_mdot, alpha, engine_onset] (auto t)
    {
        return t < engine_onset
        ? envelop_mdot * std::pow(t / unit_time(1.0), alpha)
        : engine_mdot;
    };
}

auto wind_gamma_beta(const mara::config_t& cfg, unit_scalar q)
{
    auto m0              = unit_mass(1.0);
    auto md              = unit_mass_rate(1.0);
    auto engine_onset    = unit_time(cfg.get_double("engine_onset"));
    auto engine_duration = unit_time(cfg.get_double("engine_duration"));
    auto engine_u        = unit_scalar(cfg.get_double("engine_u"));
    auto envelop_u       = unit_scalar(cfg.get_double("envelop_u"));
    auto fred            = mara::fast_rise_exponential_decay(engine_onset, engine_duration);
    auto psi             = cfg.get_double("envelop_u_index");
    auto alpha           = cfg.get_double("envelop_mdot_index");

    auto integrated_envelop_mdot = [m0, md, alpha] (unit_time t)
    {
        return m0 + md * t * std::pow(t / unit_time(1.0), alpha) / (1 + alpha);
    };

    return [integrated_envelop_mdot, envelop_u, engine_u, fred, psi, q] (auto t)
    {
        auto m0   = integrated_envelop_mdot(1.0);
        auto m    = integrated_envelop_mdot(t);
        auto p    = unit_scalar(M_PI) - q;
        auto f    = std::exp(-q * q / 0.04) + std::exp(-p * p / 0.04);
        return envelop_u * std::pow(m / m0, -psi) + engine_u * fred(t) * f;
    };
}

auto wind_profile(const mara::config_t& cfg, unit_length r, unit_scalar q, unit_time t)
{
    return mara::cold_relativistic_wind_t()
    .with_mass_loss_rate(wind_mass_loss_rate(cfg, q))
    .with_gamma_beta    (wind_gamma_beta    (cfg, q))
    .primitive(r, t);
}




//=============================================================================
solution_t initial_solution(const mara::config_t& cfg)
{
    auto r0 = unit_length(1.0);
    auto r1 = r0 * cfg.get_double("router");
    auto p0 = [cfg] (unit_length r, unit_scalar q) { return wind_profile(cfg, r, q, 1.0); };
    auto nt = cfg.get_int("nt");

    auto tween = [f=cfg.get_double("focus")] (auto x)
    {
        return M_PI * (x + f * (x - 1) * x) / (1 + 2 * f * (x - 1) * x);
    };

    auto tr = nd::linspace(0.0, 1.0, nt + 1)
    | nd::map(tween)
    | nd::adjacent_zip()
    | nd::map(util::apply_to(std::bind(sedov::generate_radial_track, r0, r1, _1, _2)));
    auto u0 = tr | nd::map(std::bind(sedov::generate_conserved, _1, p0));

    return {
        0,
        1.0,
        nd::to_shared(tr),
        nd::to_shared(u0),
    };
}

sedov::radial_godunov_data_t inner_bc(
    sedov::radial_track_t track,
    sedov::primitive_function_t primitive_func,
    nd::shared_array<srhd::primitive_t, 1> pc)
{
    auto vel  = std::min(srhd::velocity_1(front(pc)), unit_velocity(0.0));
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_moving_face_t{vel};
    auto nhat = geometric::unit_vector_on(1);
    auto pl = primitive_func(front(track.face_radii), cell_center_theta(track));
    auto pr = front(pc);
    return std::pair(srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode), vel);
}

sedov::radial_godunov_data_t outer_bc(
    sedov::radial_track_t track,
    sedov::primitive_function_t primitive_func,
    nd::shared_array<srhd::primitive_t, 1> pc)
{
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_across_contact_t();
    auto nhat = geometric::unit_vector_on(1);
    auto pl = back(pc);
    auto pr = back(pc);
    return srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode);
}

solution_t remesh(solution_t solution, unit_scalar aspect)
{
    auto t0 = solution.tracks;
    auto uc = solution.conserved;

    auto [t1, u1] = nd::unzip(nd::zip(t0, uc)
    | nd::map(util::apply_to(std::bind(sedov::refine, _1, _2, aspect))));

    return {
        solution.iteration,
        solution.time,
        t1 | nd::to_shared(),
        u1 | nd::to_shared()};
}




//=============================================================================
solution_t advance(const mara::config_t& cfg, solution_t solution, unit_time dt)
{
    auto wind = std::bind(wind_profile, cfg, _1, _2, solution.time);
    auto ibc = std::bind(inner_bc, _1, sedov::primitive_function_t{wind}, _2);
    auto obc = std::bind(outer_bc, _1, sedov::primitive_function_t{wind}, _2);

    auto t0 = solution.tracks;
    auto uc = solution.conserved;
    auto pc = nd::zip(t0, uc)             | nd::map(util::apply_to(sedov::recover_primitive));
    auto fi = nd::zip(t0, pc)             | nd::map(util::apply_to(ibc));
    auto fo = nd::zip(t0, pc)             | nd::map(util::apply_to(obc));
    auto dc = nd::zip(t0, pc)             | nd::map(util::apply_to(sedov::radial_gradient));
    auto ff = nd::zip(t0, pc, dc, fi, fo) | nd::map(util::apply_to(sedov::radial_godunov_data));
    auto dr = ff                          | nd::map(std::bind(sedov::delta_face_positions, _1, dt));

    auto gf = nd::zip(t0, pc, dc)
    | nd::extend_uniform(sedov::track_data_t())
    | nd::adjacent_zip4()
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

    return {
        solution.iteration + 1,
        solution.time + dt,
        t1 | nd::to_shared(),
        u1 | nd::to_shared(),
    };
}




//=============================================================================
auto quad_mesh(nd::shared_array<sedov::radial_track_t, 1> tracks)
{
    using vec3d = geometric::euclidean_vector_t<unit_length>;

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
void output_vtk(solution_t solution, unsigned count)
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
    auto tfinal = unit_time  (cfg.get_double("tfinal"));
    auto cfl    = unit_scalar(cfg.get_double("cfl"));
    auto aspect = cfg.get_double("aspect");
    auto vtk_it = cfg.get_int("vtk");
    auto rk     = cfg.get_int("rk");

    auto vtk_count = 0;
    output_vtk(solution, vtk_count);

    while (solution.time < tfinal)
    {
        auto dt         = cfl * nd::min(solution.tracks | nd::map(sedov::minimum_spacing)) / srhd::light_speed;
        auto num_cells  = nd::sum(solution.tracks | nd::map([] (auto t) { return size(t.face_radii); }));
        auto start      = std::chrono::high_resolution_clock::now();
        auto advance_rk = control::advance_runge_kutta(std::bind(advance, cfg, _1, dt), rk);

        solution = remesh(advance_rk(solution), aspect);

        auto stop = std::chrono::high_resolution_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        auto kzps = double(num_cells) / ms;

        std::printf("[%06lu] t=%.4f zones: %05lu kzps=%.02lf\n",
            long(solution.iteration),
            solution.time.value,
            num_cells,
            kzps);

        if (long(solution.iteration) % vtk_it == 0)
        {
            output_vtk(solution, vtk_count++);            
        }
    }
    return 0;
}
