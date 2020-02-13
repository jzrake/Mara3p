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




#pragma once
#include <cmath>
#include "core_dimensional.hpp"
#include "core_dimensional_math.hpp"
#include "core_geometric.hpp"
#include "core_numeric_tuple.hpp"




//=============================================================================
namespace iso2d {


//=============================================================================
using namespace dimensional;
using primitive_t = numeric::tuple_t<unit_mass_per_area, unit_velocity, unit_velocity>;
using conserved_t = numeric::tuple_t<unit_mass, unit_momentum, unit_momentum>;
using conserved_density_t = decltype(conserved_t() / unit_area());
using flux_vector_t       = decltype(conserved_density_t() * unit_velocity());




//=============================================================================
inline primitive_t primitive(unit_mass_per_area d, unit_velocity u, unit_velocity v)
{
    return {{d, u, v}};
}




//=============================================================================
inline auto mass_density(primitive_t p) { return numeric::get<0>(p); }
inline auto velocity_1  (primitive_t p) { return numeric::get<1>(p); }
inline auto velocity_2  (primitive_t p) { return numeric::get<2>(p); }


//=============================================================================
inline auto conserved_mass      (conserved_t u) { return numeric::get<0>(u); }
inline auto conserved_momentum_1(conserved_t u) { return numeric::get<1>(u); }
inline auto conserved_momentum_2(conserved_t u) { return numeric::get<2>(u); }


//=============================================================================
inline auto conserved_mass_density      (conserved_density_t u) { return numeric::get<0>(u); }
inline auto conserved_momentum_density_1(conserved_density_t u) { return numeric::get<1>(u); }
inline auto conserved_momentum_density_2(conserved_density_t u) { return numeric::get<2>(u); }


//=============================================================================
inline auto velocity_vector(primitive_t p)
{
    return geometric::euclidean_vector(
        velocity_1(p),
        velocity_2(p),
        unit_velocity(0.0));
}


//=============================================================================
inline auto gas_pressure(primitive_t p, unit_specific_energy cs2)
{
    return mass_density(p) * cs2;
}


//=============================================================================
inline auto outer_wavespeeds(primitive_t p, geometric::unit_vector_t nhat, unit_specific_energy cs2)
{
    auto vn = dot(nhat, velocity_vector(p));
    auto cs = sqrt(cs2);
    return std::pair(vn - cs, vn + cs);
}


//=============================================================================
inline primitive_t recover_primitive(conserved_density_t u)
{
    return numeric::tuple(
        conserved_mass_density(u),
        conserved_momentum_density_1(u) / conserved_mass_density(u),
        conserved_momentum_density_2(u) / conserved_mass_density(u));
}


//=============================================================================
inline conserved_density_t conserved_density(primitive_t p)
{
    return numeric::tuple(
        mass_density(p),
        mass_density(p) * velocity_1(p),
        mass_density(p) * velocity_2(p));
}


//=============================================================================
inline flux_vector_t flux(primitive_t p, geometric::unit_vector_t nhat, unit_specific_energy cs2)
{
    auto pg = gas_pressure(p, cs2);
    auto vn = dot(nhat, velocity_vector(p));
    auto md = unit_mass_per_area(0.0) * unit_velocity(0.0);

    auto pressure_term = numeric::tuple(
        md,
        pg * nhat.component_1(),
        pg * nhat.component_2());

    return vn * conserved_density(p) + pressure_term;
}


//=============================================================================
inline flux_vector_t riemann_hlle(
    primitive_t pl,
    primitive_t pr,
    unit_specific_energy cs2,
    geometric::unit_vector_t face_normal)
{
    auto ul = conserved_density(pl);
    auto ur = conserved_density(pr);
    auto fl = flux(pl, face_normal, cs2);
    auto fr = flux(pr, face_normal, cs2);

    auto [alm, alp] = outer_wavespeeds(pl, face_normal, cs2);
    auto [arm, arp] = outer_wavespeeds(pr, face_normal, cs2);
    auto ap = std::max(unit_velocity(0), std::max(alp, arp));
    auto am = std::min(unit_velocity(0), std::min(alm, arm));

    return (ap * fl - am * fr - (ul - ur) * ap * am) / (ap - am);
}

} // namespace iso2d




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_iso2d()
{
}

#endif // DO_UNIT_TESTS
