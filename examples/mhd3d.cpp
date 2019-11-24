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
static const double gamma_law_index = 5. / 3;
using position_t = geometric::euclidean_vector_t<dimensional::unit_length>;




//=============================================================================
struct solution_t
{
    dimensional::unit_time                       time;
    rational::number_t                           iteration;
    nd::shared_array<mhd::conserved_t,        3> conserved;
    nd::shared_array<mhd::unit_magnetic_flux, 3> magnetic_flux_1;
    nd::shared_array<mhd::unit_magnetic_flux, 3> magnetic_flux_2;
    nd::shared_array<mhd::unit_magnetic_flux, 3> magnetic_flux_3;
};




//=============================================================================
namespace h5 {

void write(const Group& group, std::string name, const solution_t& solution)
{
    write(group, "time", solution.time);
    write(group, "iteration", solution.iteration);
    write(group, "conserved", solution.conserved);
    write(group, "magnetic_flux_1", solution.magnetic_flux_1);
    write(group, "magnetic_flux_2", solution.magnetic_flux_2);
    write(group, "magnetic_flux_3", solution.magnetic_flux_3);
}

}




//=============================================================================
mhd::primitive_t initial_condition(position_t x, mhd::magnetic_field_vector_t b)
{
    return mhd::primitive(1.0, {}, 1.0, b);
}

mhd::vector_potential_t initial_vector_potential(position_t x)
{
    return mhd::vector_potential_t{};
}




//=============================================================================
solution_t initial_solution()
{
    using namespace std::placeholders;
    using geometric::component;
    using util::compose;
    using util::apply_to;

    auto N = 32;

    auto dl = dimensional::unit_length(1.0);
    auto da = dimensional::unit_area(1.0);
    auto dv = pow<3>(dimensional::unit_length(1.0) / double(N));

    auto A = [] (unsigned axis) { return compose(component(axis), initial_vector_potential); };
    auto p2c = std::bind(mhd::conserved_density, _1, gamma_law_index);

    auto xv = mesh::unit_lattice<dimensional::unit_length>(N + 1, N + 1, N + 1);
    auto xc = mesh::cell_positions(xv);
    auto [xe1, xe2, xe3] = mesh::edge_positions(xv);
    auto [ae1, ae2, ae3] = std::tuple(xe1 | nd::map(A(1)), xe2 | nd::map(A(2)), xe3 | nd::map(A(3)));
    auto [bf1, bf2, bf3] = mesh::solenoidal_difference(ae1 * dl, ae2 * dl, ae3 * dl);
    auto bc = mesh::face_to_cell(bf1, bf2, bf3) / da;
    auto pc = nd::zip(xc, bc) | nd::map(apply_to(initial_condition));
    auto uc = pc | nd::map(p2c) | nd::multiply(dv);

    return {
        dimensional::unit_time(0.0),
        rational::number(0),
        uc  | nd::to_shared(),
        bf1 | nd::to_shared(),
        bf2 | nd::to_shared(),
        bf3 | nd::to_shared(),
    };
}




//=============================================================================
int main()
{
    return 0;
}
