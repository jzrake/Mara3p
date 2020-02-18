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




#include "app_config.hpp"
#include "app_control.hpp"
#include "app_filesystem.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_config.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_rational.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"
#include "physics_iso2d.hpp"
#include "physics_two_body.hpp"
#include "scheme_plm_gradient.hpp"




//=============================================================================
auto config_template()
{
    return mara::config_template()
    .item("n",                    300)   // number of zones per side
    .item("domain_radius",        1.5)   // half-size of the domain
    .item("softening_length",    0.01)   // gravitational softening length
    .item("sink_radius",         0.01)   // radius of mass sinks
    .item("sink_rate",            1e3)   // rate of mass and momentum removal at sinks
    .item("buffer_damping_rate",  1e2)   // maximum rate of buffer driving
    .item("tfinal",             100.0)   // time to stop the simulation
    .item("restart", std::string(""))    // a checkpoint file to restart from
    .item("outdir",  std::string("."));  // the directory where output files are written
}




//=============================================================================
using namespace dimensional;
static auto mach_number = 10.0;
static auto omega_frame = dimensional::quantity_t<0, 0, -1>(1.0);




//=============================================================================
struct solution_t
{
    rational::number_t                              iteration;
    dimensional::unit_time                          time;
    nd::shared_array<iso2d::conserved_density_t, 2> conserved;
};




//=============================================================================
namespace h5 {

void write(const Group& group, std::string name, const solution_t& solution)
{
    auto sgroup = group.require_group(name);
    write(sgroup, "iteration", solution.iteration);
    write(sgroup, "time", solution.time);
    write(sgroup, "conserved", solution.conserved);
}

}




//=============================================================================
auto initial_primitive(const mara::config_t& cfg)
{
    auto softening_length = unit_length(cfg.get_double("softening_length"));

    return [softening_length] (numeric::array_t<unit_length, 2> p)
    {
        auto [x, y] = as_tuple(p);
        auto ph = two_body::potential(two_body::point_mass_t(), x, y, softening_length);
        auto r0 = sqrt(x * x + y * y);
        auto vp = sqrt(-ph) - omega_frame * r0;
        auto vx = vp * (-y / r0);
        auto vy = vp * ( x / r0);
        return iso2d::primitive(1.0, vx, vy);
    };
}

auto sound_speed_squared(numeric::array_t<unit_length, 2> p, unit_length softening_length)
{
    auto [x, y] = as_tuple(p);
    auto ph = two_body::potential(two_body::point_mass_t(), x, y, softening_length);
    return -ph / mach_number / mach_number;
}

auto buffer_rate_field(const mara::config_t& cfg)
{
    auto domain_radius = dimensional::quantity_t<0, 1, 0>(cfg.get_double("domain_radius"));
    auto buffer_rate   = dimensional::quantity_t<0, 0,-1>(cfg.get_double("buffer_damping_rate"));
    auto tightness     = dimensional::quantity_t<0,-1, 0>(5.0);

    return [domain_radius, buffer_rate, tightness] (numeric::array_t<unit_length, 2> p)
    {
        auto r = sqrt(sum(p * p));
        auto y = tightness * (r - domain_radius);
        return 0.5 * buffer_rate * (1.0 + std::tanh(y));
    };
}

auto sink_rate_field(const mara::config_t& cfg, numeric::array_t<unit_length, 2> sink_position)
{
    auto sink_radius = dimensional::quantity_t<0, 1, 0>(cfg.get_double("sink_radius"));
    auto sink_rate   = dimensional::quantity_t<0, 0,-1>(cfg.get_double("sink_rate"));

    return [sink_radius, sink_rate, sink_position] (numeric::array_t<unit_length, 2> p)
    {
        auto r6 = pow<3>(sum((p - sink_position) * (p - sink_position)));
        auto s6 = pow<6>(sink_radius);
        return sink_rate * std::exp(-r6 / s6);
    };
}

auto coriolis_term(geometric::euclidean_vector_t<unit_velocity> v)
{
    auto [vx, vy, vz] = as_tuple(v);
    return 2.0 * omega_frame * numeric::array(vy, -vx);
}

auto centrifugal_term(numeric::array_t<unit_length, 2> p)
{
    auto [x, y] = as_tuple(p);
    return omega_frame * omega_frame * numeric::array(x, y);
}


//=============================================================================
auto face_coordinates(uint n, unit_length domain_radius, uint axis)
{
    auto xv = nd::linspace(-domain_radius, domain_radius, n + 1);
    auto yv = nd::linspace(-domain_radius, domain_radius, n + 1);
    auto xc = xv | nd::adjacent_mean(0);
    auto yc = xv | nd::adjacent_mean(0);
    auto vec2 = nd::map(util::apply_to([] (auto x, auto y) { return numeric::array(x, y); }));

    switch (axis)
    {
        case 1: return nd::cartesian_product(xv, yc) | vec2 | nd::to_shared();
        case 2: return nd::cartesian_product(xc, yv) | vec2 | nd::to_shared();
    }
    throw std::invalid_argument("face_coordinates (invalid axis)");
}

auto cell_coordinates(uint n, unit_length domain_radius)
{
    auto xv = nd::linspace(-domain_radius, domain_radius, n + 1) | nd::adjacent_mean();
    auto yv = nd::linspace(-domain_radius, domain_radius, n + 1) | nd::adjacent_mean();
    auto vec2 = nd::map(util::apply_to([] (auto x, auto y) { return numeric::array(x, y); }));

    return nd::cartesian_product(xv, yv) | vec2;
}




//=============================================================================
solution_t initial(const mara::config_t& cfg)
{
    auto n = cfg.get_int("n");
    auto domain_radius = unit_length(cfg.get_double("domain_radius"));

    auto xc = cell_coordinates(n, domain_radius);
    auto u = xc | nd::map(initial_primitive(cfg)) | nd::map(iso2d::conserved_density);
    return solution_t{0, 0.0, nd::to_shared(u)};
}

solution_t next(const mara::config_t& cfg, solution_t solution)
{
    using namespace std::placeholders;

    auto n                = cfg.get_int("n");
    auto domain_radius    = unit_length(cfg.get_double("domain_radius"));
    auto softening_length = unit_length(cfg.get_double("softening_length"));

    auto riemann = [rs=softening_length] (uint axis)
    {
        return [axis, rs] (auto pl, auto pr, auto xf)
        {
            auto cs2 = sound_speed_squared(xf, rs);
            return iso2d::riemann_hlle(pl, pr, cs2, geometric::unit_vector_on(axis));
        };
    };

    auto elements = two_body::orbital_elements();
    auto inertial = two_body::orbital_state(elements, solution.time);
    auto state    = two_body::rotate(inertial, omega_frame * solution.time);

    auto gravitational_acceleration = [softening_length] (auto component)
    {
        return [c=component, s=softening_length] (auto p)
        {
            return two_body::gravitational_acceleration(c, p[0], p[1], s);
        };
    };

    auto xf1 = face_coordinates(n, domain_radius, 1);
    auto xf2 = face_coordinates(n, domain_radius, 2);
    auto xc  = cell_coordinates(n, domain_radius);
    auto u0  = solution.conserved;
    auto pc  = u0 | nd::map(iso2d::recover_primitive) | nd::to_shared();

    auto gx = pc | nd::extend_zero_gradient(0) | nd::adjacent_zip3(0) | nd::map(mara::plm_gradient(2.0)) | nd::to_shared();
    auto gy = pc | nd::extend_zero_gradient(1) | nd::adjacent_zip3(1) | nd::map(mara::plm_gradient(2.0)) | nd::to_shared();
    auto plx = pc - 0.5 * gx;
    auto prx = pc + 0.5 * gx;
    auto ply = pc - 0.5 * gy;
    auto pry = pc + 0.5 * gy;

    auto fx = nd::zip(
        prx | nd::select(0, 0, -1),
        plx | nd::select(0, 1),
        xf1 | nd::select(0, 1, -1))
    | nd::map(util::apply_to(riemann(1)))
    | nd::extend_zero_gradient(0)
    | nd::to_shared();

    auto fy = nd::zip(
        pry | nd::select(1, 0, -1),
        ply | nd::select(1, 1),
        xf2 | nd::select(1, 1, -1))
    | nd::map(util::apply_to(riemann(2)))
    | nd::extend_zero_gradient(1)
    | nd::to_shared();

    auto dl  = 2.0 * domain_radius / double(n);
    auto dt  = 0.25 * dl / nd::max(pc | nd::map(std::bind(iso2d::fastest_wavespeed, _1, unit_specific_energy(0.0))));

    auto dfx = fx | nd::adjacent_diff(0);
    auto dfy = fy | nd::adjacent_diff(1);

    auto ag1 = xc | nd::map(gravitational_acceleration(state.first));
    auto ag2 = xc | nd::map(gravitational_acceleration(state.second));
    auto cen = xc | nd::map(centrifugal_term);
    auto cor = pc | nd::map(iso2d::velocity_vector) | nd::map(coriolis_term);

    auto ss1 = -u0 * (xc | nd::map(sink_rate_field(cfg, two_body::position(state.first))));
    auto ss2 = -u0 * (xc | nd::map(sink_rate_field(cfg, two_body::position(state.second))));
    auto sg1 = nd::zip(pc, ag1) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));
    auto sg2 = nd::zip(pc, ag2) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));
    auto sf1 = nd::zip(pc, cen) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));
    auto sf2 = nd::zip(pc, cor) | nd::map(util::apply_to(iso2d::acceleration_to_source_terms));

    auto uinit = xc | nd::map(initial_primitive(cfg)) | nd::map(iso2d::conserved_density);
    auto br    = xc | nd::map(buffer_rate_field(cfg));

    auto sb = -(u0 - uinit) * br;
    auto u1 = u0 - (dfx + dfy) * dt / dl + (sg1 + sg2 + ss1 + ss2 + sf1 + sf2 + sb) * dt;

    return {
        solution.iteration + 1,
        solution.time + dt,
        u1 | nd::to_shared()
    };
}




//=============================================================================
void side_effects(solution_t solution)
{
    if (long(solution.iteration) % 50 == 0)
    {
        auto fname = util::format("chkpt.%04ld.h5", long(solution.iteration) / 50);
        auto h5f = h5::File(fname, "w");
        h5::write(h5f, "solution", solution);
        std::printf("write %s\n", fname.data());
    }
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto cfg = config_template().create().update(mara::argv_to_string_map(argc, argv));
    auto solution = initial(cfg);

    auto ncells = cfg.get_int("n") * cfg.get_int("n");
    auto tfinal = unit_time(cfg.get_double("tfinal"));

    while (solution.time < tfinal)
    {
        side_effects(solution);

        auto [s1, ticks] = control::invoke_timed(next, cfg, solution);
        auto kzps = 1e6 * ncells / std::chrono::duration_cast<std::chrono::nanoseconds>(ticks).count();

        solution = s1;

        std::printf("[%ld] t=%lf kzps=%lf\n", long(solution.iteration), solution.time.value, kzps);
    }

    return 0;
}
