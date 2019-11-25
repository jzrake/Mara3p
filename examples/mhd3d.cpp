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




#include <cstdio>
#include "app_config.hpp"
#include "app_control.hpp"
#include "app_hdf5.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_geometric.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_rational.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_sequence.hpp"
#include "physics_mhd.hpp"
#include "mesh_cartesian_3d.hpp"




//=============================================================================
static const unsigned zone_count = 64;
static const double gamma_law_index = 5. / 3;




//=============================================================================
struct solution_t
{
    rational::number_t                           iteration;
    dimensional::unit_time                       time;
    nd::shared_array<mhd::conserved_t,        3> conserved;
    nd::shared_array<mhd::unit_magnetic_flux, 3> magnetic_flux_1;
    nd::shared_array<mhd::unit_magnetic_flux, 3> magnetic_flux_2;
    nd::shared_array<mhd::unit_magnetic_flux, 3> magnetic_flux_3;
};




//=============================================================================
using namespace std::placeholders;
using util::apply_to;
using util::compose;
using geometric::unit_vector_on;
using geometric::component;
using position_t            = geometric::euclidean_vector_t<dimensional::unit_length>;
using solution_with_tasks_t = std::pair<solution_t, control::task_t>;
using timed_state_pair_t    = control::timed_pair_t<solution_with_tasks_t>;




auto magnetic_flux(solution_t s) { return std::tuple(s.magnetic_flux_1, s.magnetic_flux_2, s.magnetic_flux_3); }
auto conserved(solution_t s) { return s.conserved; }
auto solution(solution_with_tasks_t p) { return p.first; }
auto tasks(solution_with_tasks_t p)   { return p.second; }

template<typename P>
void print_shape(const char* label, nd::array_t<P, 3> a)
{
    std::printf("%s: %lu %lu %lu\n", label, a.shape[0], a.shape[1], a.shape[2]);
}




//=============================================================================
namespace h5 {

void write(const Group& group, const solution_t& solution)
{
    write(group, "iteration", solution.iteration);
    write(group, "time", solution.time);
    write(group, "conserved", solution.conserved);
    write(group, "magnetic_flux_1", solution.magnetic_flux_1);
    write(group, "magnetic_flux_2", solution.magnetic_flux_2);
    write(group, "magnetic_flux_3", solution.magnetic_flux_3);
}

}




//=============================================================================
mhd::primitive_t initial_primitive(position_t x, mhd::magnetic_field_vector_t b)
{
    auto x0 = geometric::euclidean_vector<dimensional::unit_length>(0.5, 0.5, 0.5);
    auto R = sqrt(length_squared(x - x0));
    auto d = R < dimensional::unit_length(0.125) ? 1.0 : 0.1;
    auto p = R < dimensional::unit_length(0.125) ? 1.0 : 0.125;
    return mhd::primitive(d, {}, p, b);
}

mhd::vector_potential_t initial_vector_potential(position_t x)
{
    return mhd::vector_potential_t{0.0, x.component_1().value * 4.0, 0.0};
}




//=============================================================================
solution_t initial_solution()
{
    auto N = zone_count;

    auto dl = pow<1>(dimensional::unit_length(1.0) / double(N));
    auto da = pow<2>(dimensional::unit_length(1.0) / double(N));
    auto dv = pow<3>(dimensional::unit_length(1.0) / double(N));

    auto A = [] (unsigned dir) { return compose(component(dir), initial_vector_potential); };
    auto p2c = std::bind(mhd::conserved_density, _1, gamma_law_index);

    auto xv = mesh::unit_lattice<dimensional::unit_length>(N + 1, N + 1, N + 1);
    auto xc = mesh::cell_positions(xv);
    auto [xe1, xe2, xe3] = mesh::edge_positions(xv);
    auto [ae1, ae2, ae3] = std::tuple(xe1 | nd::map(A(1)), xe2 | nd::map(A(2)), xe3 | nd::map(A(3)));
    auto [mf1, mf2, mf3] = mesh::solenoidal_difference(ae1 * dl, ae2 * dl, ae3 * dl);
    auto bc = mesh::face_to_cell(mf1, mf2, mf3) / da;
    auto pc = nd::zip(xc, bc) | nd::map(apply_to(initial_primitive));
    auto uc = pc | nd::map(p2c) | nd::multiply(dv);

    return {
        rational::number(0),
        dimensional::unit_time(0.0),
        uc  | nd::to_shared(),
        mf1 | nd::to_shared(),
        mf2 | nd::to_shared(),
        mf3 | nd::to_shared(),
    };
}




// These are riemann solver functions for each of the three axes.
auto rs1 = apply_to(std::bind(mhd::riemann_hlle, _1, _2, _3, unit_vector_on(1), gamma_law_index));
auto rs2 = apply_to(std::bind(mhd::riemann_hlle, _1, _2, _3, unit_vector_on(2), gamma_law_index));
auto rs3 = apply_to(std::bind(mhd::riemann_hlle, _1, _2, _3, unit_vector_on(3), gamma_law_index));




// These are "remove-transverse" functions; they return a rectangular pipe from
// the middle of a 3D array along a given axis.
auto rt1 = mesh::remove_transverse_i(1);
auto rt2 = mesh::remove_transverse_j(1);
auto rt3 = mesh::remove_transverse_k(1);




// These are "remove-longitudinal" functions; they remove a layer of zones from
// the ends of a 3D array along a given axis.
auto rl1 = nd::select(0, 1, -1);
auto rl2 = nd::select(1, 1, -1);
auto rl3 = nd::select(2, 1, -1);




template<typename P, typename Q>
auto godunov_fluxes(nd::array_t<P, 3> primitive_array, nd::array_t<Q, 3> bf1, nd::array_t<Q, 3> bf2, nd::array_t<Q, 3> bf3)
{
    auto bf1_ = bf1 | nd::extend_periodic(1) | nd::extend_periodic(2);
    auto bf2_ = bf2 | nd::extend_periodic(2) | nd::extend_periodic(0);
    auto bf3_ = bf3 | nd::extend_periodic(0) | nd::extend_periodic(1);

    auto p0 = primitive_array | mesh::extend_periodic(1) | nd::to_shared();

    auto pf1 = nd::zip(p0 | nd::select(0, 0, -1), p0 | nd::select(0, 1), bf1_);
    auto pf2 = nd::zip(p0 | nd::select(1, 0, -1), p0 | nd::select(1, 1), bf2_);
    auto pf3 = nd::zip(p0 | nd::select(2, 0, -1), p0 | nd::select(2, 1), bf3_);

    auto [ff1, ef1] = nd::unzip(pf1 | nd::map(rs1) | nd::to_shared());
    auto [ff2, ef2] = nd::unzip(pf2 | nd::map(rs2) | nd::to_shared());
    auto [ff3, ef3] = nd::unzip(pf3 | nd::map(rs3) | nd::to_shared());

    return std::tuple(ff1, ff2, ff3, ef1, ef2, ef3);
}




template<typename P>
auto edge_emf_from_face(nd::array_t<P, 3> ef1, nd::array_t<P, 3> ef2, nd::array_t<P, 3> ef3)
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




//=============================================================================
solution_t advance(solution_t solution)
{
    auto c2p = apply_to(std::bind(mhd::recover_primitive, _1, _2, gamma_law_index));
    auto N = zone_count;

    auto ev = nd::to_shared();
    auto dl = pow<1>(dimensional::unit_length(1.0) / double(N));
    auto da = pow<2>(dimensional::unit_length(1.0) / double(N));
    auto dv = pow<3>(dimensional::unit_length(1.0) / double(N));
    auto dt = 0.02 * dl / dimensional::unit_velocity(1.0);

    auto [mf1, mf2, mf3] = magnetic_flux(solution);
    auto [bf1, bf2, bf3] = std::tuple(mf1 / da, mf2 / da, mf3 / da);
    auto uc = conserved(solution);
    auto bc = mesh::face_to_cell(mf1, mf2, mf3);
    auto pc = zip(uc / dv, bc / da) | nd::map(c2p) | nd::to_shared();

    auto [ff1, ff2, ff3, ef1, ef2, ef3] = godunov_fluxes(pc, bf1, bf2, bf3);
    auto [ee1, ee2, ee3] = edge_emf_from_face(ef1, ef2, ef3);
    auto [cf1, cf2, cf3] = mesh::solenoidal_difference(ee1|rl1|ev, ee2|rl2|ev, ee3|rl3|ev);
    auto dc              = mesh::divergence_difference(ff1|rt1|ev, ff2|rt2|ev, ff3|rt3|ev);

    auto uc_  = uc  - dc  * da * dt;
    auto mf1_ = mf1 - cf1 * dl * dt;
    auto mf2_ = mf2 - cf2 * dl * dt;
    auto mf3_ = mf3 - cf3 * dl * dt;

    return {
        solution.iteration + 1,
        solution.time + dt,
        uc_  | nd::to_shared(),
        mf1_ | nd::to_shared(),
        mf2_ | nd::to_shared(),
        mf3_ | nd::to_shared(),
    };
}

solution_with_tasks_t advance_app_state(solution_with_tasks_t state)
{
    return std::pair(advance(solution(state)), jump(tasks(state), solution(state).time, dimensional::unit_time(0.01)));
}

solution_with_tasks_t initial_app_state(const mara::config_t& run_config)
{
    return std::pair(initial_solution(), control::task("write_checkpoint", 0.0));
}

bool should_continue(timed_state_pair_t state_pair)
{
    auto this_soln = solution(control::this_state(state_pair));

    return long(this_soln.iteration) < 100;
}




//=============================================================================
void write_checkpoint(solution_t solution, long count)
{
    auto fname = util::format("mhd3d.%04lu.h5", count);
    std::printf("Write checkpoint %s\n", fname.data());
    write(h5::File(fname, "w"), solution);
}

void print_run_loop(timed_state_pair_t p)
{
    auto soln = solution(control::this_state(p));
    auto nz = size(soln.conserved);
    auto us = control::microseconds_separating(p);

    std::printf("[%07lu] t=%.3lf dt=%.2e zones=%lu Mzps=%.2lf\n",
        long(soln.iteration),
        soln.time.value,
        0.2,
        nz,
        nz / us);
}

void side_effects(timed_state_pair_t p)
{
    auto this_soln = solution(control::this_state(p));
    auto last_soln = solution(control::last_state(p));
    auto this_task = tasks(control::this_state(p));
    auto last_task = tasks(control::last_state(p));

    if (this_task.count != last_task.count)
        write_checkpoint(last_soln, last_task.count);

    print_run_loop(p);
}




//=============================================================================
int main()
{
    auto run_config = mara::config_t();

    auto simulation = seq::generate(initial_app_state(run_config), advance_app_state)
    | seq::pair_with(control::time_point_sequence())
    | seq::window()
    | seq::take_while(should_continue);

    for (auto state_pair : simulation)
    {
        side_effects(state_pair);
    }
    return 0;
}
