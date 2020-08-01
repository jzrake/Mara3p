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
#include "mara.hpp"
#include "app_config.hpp"
#include "app_control.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_config.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_rational.hpp"
#include "app_filesystem.hpp"
#include "app_vtk.hpp"
#include "core_util.hpp"
#include "model_wind.hpp"
#include "parallel_thread_pool.hpp"
#include "scheme_sedov2d.hpp"

using namespace dimensional;
using namespace std::placeholders;
static mara::ThreadPool thread_pool;




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("nt",                   128, "number of radial tracks")
    .item("rk",                     1, "Runge-Kutta order (1, 2, or 3)")
    .item("max_aspect",           2.0, "aspect ratio of the longest cell allowed")
    .item("min_aspect",           0.5, "aspect ratio of the widest cell allowed")
    .item("inner_velocity",       0.0, "the velocity with which to move the inner boundary")
    .item("outer_velocity",       0.9, "the boundary boundary velocity (also determines the outer radius)")
    .item("focus",                0.0, "amount by which to increase resolution near the poles [0, 1)")
    .item("tfinal",               1e3, "time to stop the simulation")
    .item("tstart",               2.0, "time (following the merger) when the simulation begins")
    .item("cpi",                 5000, "number of iterations between checkpoints")
    .item("vtk",                    0, "number of iterations between VTK outputs")
    .item("cfl",                  0.5, "courant number")
    .item("threads",                1, "the number of concurrent threads to execute on")
    .item("restart",               "", "a checkpoint file to restart from")
    .item("outdir",               ".", "the directory where output files are written")
    .item("envelop_mass"   , 1e-3, "mass of the relativistic envelop (solar mass)")
    .item("cloud_mass"     , 1e-2, "mass of the merger ejecta cloud (solar mass)")
    .item("psi"            , 0.22, "index psi in u(m) ~ m^-psi")
    .item("u_min"          , 0.05, "four-velocity of the slowest envelop shell")
    .item("engine_delay"   , 1e3,  "time following the cloud onset when the jet begins")
    .item("engine_duration", 1e3,  "duration of the engine")
    .item("engine_theta"   , 0.10, "engine opening angle")
    .item("engine_u"       , 10.0, "engine four-velocity")
    .item("engine_strength", 0.01, "E / M c^2: M = cloud mass, E = isotropic-equivalent jet energy");
}




//=============================================================================
template<typename P>
auto evaluate_on(mara::ThreadPool& pool, nd::array_t<P, 1> array)
{
    using value_type = std::remove_cv_t<typename nd::array_t<P, 1>::value_type>;
    auto futures = nd::make_unique_array<std::future<value_type>>(shape(array));
    auto results = nd::make_unique_array<value_type>(shape(array));

    for (std::size_t i = 0; i < size(array); ++i)
    {
        futures(i) = pool.enqueue(array, i);
    }
    for (std::size_t i = 0; i < size(array); ++i)
    {
        results(i) = futures(i).get();
    }
    return nd::make_shared_array(std::move(results));
}

auto evaluate_on(mara::ThreadPool& pool)
{
    return [&pool] (auto array)
    {
        return evaluate_on(pool, array);
    };
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
    auto eval   = evaluate_on(thread_pool);
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
    })) | eval;
    result.conserved = nd::zip(s.conserved, t.conserved)
    | nd::map(util::apply_to([b] (auto s, auto t)
    {
        return nd::to_shared(s * double(b) + t * double(1 - b));
    })) | eval;

    return result;
}




//=============================================================================
void write_checkpoint(const mara::config_t& cfg, solution_t solution)
{
    auto count     = long(solution.iteration) / cfg.get_int("cpi");
    auto outdir    = cfg.get_string("outdir");
    auto fname     = util::format("%s/chkpt.%04d.h5", outdir.data(), count);
    auto file      = h5::File(fname, "w");
    auto tracks    = h5::Group(file).require_group("tracks");

    std::printf("output %s\n", fname.data());

    h5::write(file, "git_commit", std::string(MARA_GIT_COMMIT));
    h5::write(file, "run_config", cfg);
    h5::write(file, "iteration", solution.iteration);
    h5::write(file, "time", solution.time);

    for (auto n : nd::range(size(solution.tracks)))
    {
        auto track = tracks.require_group(std::to_string(n));
        h5::write(track, "face_radii", solution.tracks(n).face_radii);
        h5::write(track, "theta0", solution.tracks(n).theta0);
        h5::write(track, "theta1", solution.tracks(n).theta1);
        h5::write(track, "conserved", solution.conserved(n));
        h5::write(track, "primitive", sedov::recover_primitive(solution.tracks(n), solution.conserved(n)));
    }
}

solution_t read_checkpoint(std::string fname)
{
    auto file         = h5::File(fname, "r");
    auto tracks_group = h5::Group(file).open_group("tracks");
    auto cfg          = mara::config_parameter_map_t();
    auto tracks       = nd::make_unique_array<sedov::radial_track_t>(nd::uivec(tracks_group.size()));
    auto conserved    = nd::make_unique_array<nd::shared_array<srhd::conserved_t, 1>>(nd::uivec(tracks_group.size()));
    auto solution     = solution_t();

    for (auto name : tracks_group)
    {
        auto n = std::stoi(name);
        h5::read(tracks_group.open_group(name), "face_radii", tracks(n).face_radii);
        h5::read(tracks_group.open_group(name), "theta0", tracks(n).theta0);
        h5::read(tracks_group.open_group(name), "theta1", tracks(n).theta1);
        h5::read(tracks_group.open_group(name), "conserved", conserved(n));
    }

    h5::read(file, "iteration", solution.iteration);
    h5::read(file, "time", solution.time);
    solution.tracks    = nd::make_shared_array(std::move(tracks));
    solution.conserved = nd::make_shared_array(std::move(conserved));

    return solution;
}

mara::config_parameter_map_t restart_run_config(const mara::config_string_map_t& args)
{
    if (args.count("restart"))
    {
        auto file = h5::File(args.at("restart"), "r");
        return h5::read<mara::config_parameter_map_t>(file, "run_config");
    }
    return {};
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
void write_vtk(const mara::config_t& cfg, solution_t solution)
{
    auto count     = long(solution.iteration) / cfg.get_int("vtk");
    auto outdir    = cfg.get_string("outdir");
    auto fname     = util::format("%s/primitive.%04d.vtk", outdir.data(), count);
    auto outf = std::ofstream(fname, std::ios_base::out);
    auto [vertices, indexes] = quad_mesh(solution.tracks);

    auto p0 = nd::zip(solution.tracks, solution.conserved)
    | nd::map(util::apply_to(sedov::recover_primitive))
    | nd::to_shared()
    | nd::flat();

    auto density  = p0 | nd::map(srhd::mass_density) | nd::to_shared();
    auto pressure = p0 | nd::map(srhd::gas_pressure) | nd::to_shared();
    auto ur       = p0 | nd::map(srhd::gamma_beta_1) | nd::to_shared();

    std::printf("output %s\n", fname.data());

    vtk::write(outf, "Grid", vertices, indexes,
        std::pair("density", density),
        std::pair("pressure", pressure),
        std::pair("radial-gamma-beta", ur));
}




//=============================================================================
class jet_wind_nozzle_prescription_t
{
public:


    jet_wind_nozzle_prescription_t(const mara::config_t& cfg)
    {
        envelop_mass    = cfg.get_double("envelop_mass");  
        cloud_mass      = cfg.get_double("cloud_mass");
        psi             = cfg.get_double("psi");
        u_min           = cfg.get_double("u_min");
        engine_delay    = cfg.get_double("engine_delay");
        engine_duration = cfg.get_double("engine_duration");
        engine_theta    = cfg.get_double("engine_theta");
        engine_u        = cfg.get_double("engine_u");
        engine_strength = cfg.get_double("engine_strength");
    }


    void print(std::ostream& out) const
    {
        out << "nozzle description:\n\n";
        out << "\tt1 = " << get_t1().value << '\n';
        out << "\tt2 = " << get_t2().value << '\n';
        out << "\tt3 = " << get_t3().value << '\n';
        out << std::endl;
    }


    /**
     * @brief      The time when the slowest envelop shell comes through r =
     *             1.0.
     *
     * @return     A time
     */
    unit_time get_t1() const
    {
        auto v_min = c * u_min / sqrt(1.0 + u_min * u_min);
        return r_in / v_min;
    }


    /**
     * @brief      The time when the jet turns on.
     *
     * @return     A time
     */
    unit_time get_t2() const
    {
        return get_t1() + engine_delay;
    }


    /**
     * @brief      The time when the jet turns off.
     *
     * @return     A time
     */
    unit_time get_t3() const
    {
        return get_t2() + engine_duration;
    }


    /**
     * @brief      Return the mass flux, per sterian, in the jet
     *
     * @return     The mass flux
     */
    unit_mass_rate get_engine_mass_flux() const
    {
        auto engine_gamma = sqrt(1.0 + engine_u * engine_u);
        auto E = engine_strength * cloud_mass;
        auto L = E / (4.0 * M_PI * engine_duration);
        auto Mdot = L / engine_gamma;
        return Mdot;
    }


    /**
     * @brief      Return true of a polar angle is within theta_jet of either
     *             pole
     *
     * @param[in]  q     The polar angle theta
     *
     * @return     True or false
     */
    bool in_nozzle(unit_scalar q) const
    {
        return q < engine_theta || q > unit_scalar(M_PI) - engine_theta;
    }


    /**
     * @brief      Determine the zone of the ambient medium for a given radius
     *             and time.
     *
     * @param[in]  r     The radius
     * @param[in]  t     The time
     *
     * @return     An integer
     */
    int get_zone(unit_length r, unit_time t) const
    {
        auto v_min = c * u_min / sqrt(1 + u_min * u_min);

        if (t <= r / v_min)
        {
            return 1;
        }
        else if (t <= r / v_min + engine_delay)
        {
            return 2;
        }
        else if (t <= r / v_min + engine_delay + engine_duration)
        {
            return 3;
        }
        else
        {
            return 4;            
        }
    }


    /**
     * @brief      Return the radial four-velocity (gamma-beta).
     *
     * @param[in]  r     The radius
     * @param[in]  q     The polar angle theta
     * @param[in]  t     The time
     *
     * @return     Gamma-beta
     */
    unit_scalar gamma_beta(unit_length r, unit_scalar q, unit_time t) const
    {
        if (get_zone(r, t) == 1)
        {
            auto v = r / t;
            auto b = v / c;
            auto u = b / sqrt(1 - b * b);
            return u;
        }
        else if (get_zone(r, t) == 2)
        {
            return u_min;
        }
        else if (get_zone(r, t) == 3 && in_nozzle(q))
        {
            return engine_u;
        }
        else
        {
            return u_min;
        }
    }


    /**
     * @brief      Return the mass flux per solid angle.
     *
     * @param[in]  r     The radius
     * @param[in]  q     The polar angle theta
     * @param[in]  t     The time
     *
     * @return     The mass flux per steradian
     */
    unit_mass_rate mass_flux(unit_length r, unit_scalar q, unit_time t) const
    {
        auto m0 = envelop_mass;
        auto mc = cloud_mass;

        if (get_zone(r, t) == 1)
        {
            auto s = r / c / t;
            auto f = std::pow(s / u_min, -1.0 / psi) * std::pow(1 - s * s, 0.5 / psi - 1.0);
            return m0 / (4.0 * M_PI * psi * t) * f;
        }
        else if (get_zone(r, t) == 2)
        {
            return mc / (engine_delay * 4.0 * M_PI);
        }
        else if (get_zone(r, t) == 3 && in_nozzle(q))
        {
            return get_engine_mass_flux();
        }
        else
        {
            return mc / (engine_delay * 4.0 * M_PI);
        }
    }


    /**
     * @brief      Return the primitive variable state, valid in zones 1 or 2 or
     *             at the inner boundary.
     *
     * @param[in]  r     The radius
     * @param[in]  q     The polar angle theta
     * @param[in]  t     The time
     *
     * @return     The primitive variable state
     */
    srhd::primitive_t primitive_state(unit_length r, unit_scalar q, unit_time t) const
    {
        auto mdot = mass_flux(r, q, t);
        auto u = gamma_beta(r, q, t);
        auto d = mdot / (r * r * u * c);
        auto s = 1e-6; // p / rho^(4/3)
        auto p = s * std::pow(d.value, 4. / 3);
        return srhd::primitive(d, u, 0.0, 0.0, p);
    }


private:
    unit_velocity  c               = srhd::light_speed;
    unit_length    r_in            = 1.0;  // inner boundary radius
    unit_mass      envelop_mass    = 1e-3; // mass of the relativistic envelop (solar mass)
    unit_mass      cloud_mass      = 1e-2; // mass of the merger ejecta cloud (solar mass)
    unit_scalar    psi             = 0.22; // index psi in u(m) ~ m^-psi
    unit_scalar    u_min           = 0.05; // four-velocity of the slowest envelop shell
    unit_time      engine_delay    = 1e3;  // time following the cloud onset when the jet begins
    unit_time      engine_duration = 1e3;  // duration of the engine
    unit_scalar    engine_theta    = 0.10; // engine opening angle
    unit_scalar    engine_u        = 10.0; // engine four-velocity
    unit_scalar    engine_strength = 0.01; // E / M c^2: M = cloud mass, E = isotropic-equivalent jet energy
};




//=============================================================================
auto wind_profile(const mara::config_t& cfg, unit_length r, unit_scalar q, unit_time t)
{
    return jet_wind_nozzle_prescription_t(cfg).primitive_state(r, q, t);
}




//=============================================================================
solution_t initial_solution(const mara::config_t& cfg)
{
    if (! cfg.get_string("restart").empty())
    {
        return read_checkpoint(cfg.get_string("restart"));
    }

    auto eval = evaluate_on(thread_pool);
    auto ts = unit_time(cfg.get_double("tstart"));
    auto vo = unit_velocity(cfg.get_double("outer_velocity"));
    auto r0 = unit_length(1.0);
    auto r1 = ts * vo;
    auto p0 = [cfg, ts] (unit_length r, unit_scalar q) { return wind_profile(cfg, r, q, ts); };
    auto nt = cfg.get_int("nt");
    auto aspect = cfg.get_double("max_aspect");

    auto tween = [f=cfg.get_double("focus")] (auto x)
    {
        return M_PI * (x + f * (x - 1) * x) / (1 + 2 * f * (x - 1) * x);
    };

    auto tr = nd::linspace(0.0, 1.0, nt + 1)
    | nd::map(tween)
    | nd::adjacent_zip()
    | nd::mapv(std::bind(sedov::generate_radial_track, r0, r1, _1, _2, aspect));
    auto u0 = tr | nd::map(std::bind(sedov::generate_conserved, _1, p0));

    return {
        0,
        ts,
        eval(tr),
        eval(u0),
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
    auto vel  = srhd::velocity_1(primitive_func(back(track.face_radii), cell_center_theta(track)));
    auto mode = srhd::riemann_solver_mode_hllc_fluxes_moving_face_t{vel};
    auto nhat = geometric::unit_vector_on(1);
    auto pl = back(pc);
    auto pr = primitive_func(back(track.face_radii), cell_center_theta(track));
    return std::pair(srhd::riemann_solver(pl, pr, nhat, 4. / 3, mode), vel);
}

solution_t remesh(solution_t solution, unit_scalar max_aspect, unit_scalar min_aspect, unit_length inner_boundary_radius)
{
    auto eval = evaluate_on(thread_pool);
    auto t0 = solution.tracks;
    auto uc = solution.conserved;
    auto [t1, u1] = nd::unzip(nd::zip(t0, uc) | nd::mapv(std::bind(sedov::remesh, _1, _2, max_aspect, min_aspect, inner_boundary_radius)));

    return {
        solution.iteration,
        solution.time,
        t1 | eval,
        u1 | eval,
    };
}




//=============================================================================
solution_t advance(const mara::config_t& cfg, solution_t solution, unit_time dt, bool use_plm)
{
    auto eval = evaluate_on(thread_pool);
    auto wind = std::bind(wind_profile, cfg, _1, _2, solution.time);
    auto ibc = std::bind(inner_bc, _1, sedov::primitive_function_t{wind}, _2);
    auto obc = std::bind(outer_bc, _1, sedov::primitive_function_t{wind}, _2);
    auto radial_gradient = std::bind(sedov::radial_gradient, _1, _2, use_plm);
    auto polar_godunov   = std::bind(sedov::polar_godunov_data, _1, _2, _3, _4, use_plm);

    auto t0 = solution.tracks;
    auto uc = solution.conserved;
    auto pc = nd::zip(t0, uc)             | nd::mapv(sedov::recover_primitive)       | eval;
    auto fi = nd::zip(t0, pc)             | nd::mapv(ibc)                            | nd::to_shared();
    auto fo = nd::zip(t0, pc)             | nd::mapv(obc)                            | nd::to_shared();
    auto dc = nd::zip(t0, pc)             | nd::mapv(radial_gradient)                | eval;
    auto ff = nd::zip(t0, pc, dc, fi, fo) | nd::mapv(sedov::radial_godunov_data)     | eval;
    auto dr = ff                          | nd::map(std::bind(sedov::delta_face_positions, _1, dt)) | eval;

    auto gf = nd::zip(t0, pc, dc)
    | nd::extend_uniform(sedov::track_data_t())
    | nd::adjacent_zip4()
    | nd::map(sedov::copy_track_data4)
    | nd::mapv(polar_godunov)
    | nd::extend_uniform(nd::shared_array<sedov::polar_godunov_data_t, 1>{})
    | eval;

    auto u1 = nd::range(size(uc))
    | nd::map([t0, uc, pc, ff, gf, dt] (nd::uint j)
    {
        return nd::to_shared(uc(j) + sedov::delta_conserved(t0(j), pc(j), ff(j), gf(j), gf(j + 1), dt));
    })
    | eval;

    auto t1 = nd::zip(t0, dr)
    | nd::map(util::apply_to([] (auto track, auto dr)
    {
        return sedov::radial_track_t{
            nd::to_shared(track.face_radii + dr),
            track.theta0,
            track.theta1,
        };
    }))
    | eval;

    return {
        solution.iteration + 1,
        solution.time + dt,
        t1,
        u1,
    };
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto args = mara::argv_to_string_map(argc, argv);

    auto cfg = config_template()
    .create()
    .update(restart_run_config(args))
    .update(args);

    thread_pool.reset(cfg.get_int("threads"));

    auto tfinal       = unit_time  (cfg.get_double("tfinal"));
    // auto engine_onset = unit_time  (cfg.get_double("engine_onset"));
    auto cfl          = unit_scalar(cfg.get_double("cfl"));
    // auto inner_vel    = unit_velocity(cfg.get_double("inner_velocity"));
    auto max_aspect   = cfg.get_double("max_aspect");
    auto min_aspect   = cfg.get_double("min_aspect");
    auto vtk          = cfg.get_int("vtk");
    auto cpi          = cfg.get_int("cpi");
    auto rk           = cfg.get_int("rk");
    auto outdir       = cfg.get_string("outdir");
    auto solution     = initial_solution(cfg);
    auto backup       = solution_t{};
    auto safety       = false;
    auto failed_iter  = rational::number_t();

    mara::pretty_print(std::cout, "config", cfg);
    mara::filesystem::require_dir(outdir);
    jet_wind_nozzle_prescription_t(cfg).print(std::cout);

    if (vtk > 0) write_vtk(cfg, solution);
    if (cpi > 0) write_checkpoint(cfg, solution);        

    while (solution.time < tfinal)
    {
        auto cfl_factor = 1.0;//solution.time < unit_time(cfg.get_double("engine_onset")) ? 1.0 : 4.0;
        auto dt         = cfl * (safety ? 1e-3 : 1.0) * nd::min(solution.tracks | nd::map(sedov::minimum_spacing)) / srhd::light_speed * cfl_factor;
        auto num_cells  = nd::sum(solution.tracks | nd::map([] (auto t) { return size(t.face_radii); }));
        auto start      = std::chrono::high_resolution_clock::now();
        auto advance_rk = control::advance_runge_kutta(std::bind(advance, cfg, _1, dt, ! safety), rk);
        auto inner_boundary_radius = unit_length(1.0); // std::max(unit_time(1.0), solution.time - engine_onset) * inner_vel;

        try {
            auto pre_step = solution;

            solution = remesh(advance_rk(solution), max_aspect, min_aspect, inner_boundary_radius);
            backup   = pre_step;
            safety   = false;

            auto stop = std::chrono::high_resolution_clock::now();
            auto ms   = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
            auto kzps = double(num_cells) / ms;

            std::printf("[%06lu] t=%.4f dt=%.04e zones: %05lu kzps=%.02lf\n",
                long(solution.iteration),
                solution.time.value,
                dt.value,
                num_cells,
                kzps);
            std::fflush(stdout);
        }
        catch (const std::exception& e)
        {
            if (solution.iteration == failed_iter)
            {
                throw;
            }

            failed_iter = solution.iteration;
            solution    = backup;
            safety      = true;

            std::printf("%s\n[re-trying time step %lu -> %lu in safety mode]\n",
                e.what(),
                long(solution.iteration),
                long(solution.iteration + 1));
        }

        if (vtk > 0 && long(solution.iteration) % vtk == 0)
        {
            write_vtk(cfg, solution);
        }
        if (cpi > 0 && long(solution.iteration) % cpi == 0)
        {
            write_checkpoint(cfg, solution);
        }
    }
    return 0;
}
