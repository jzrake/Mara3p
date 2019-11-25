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




#pragma once
#include <cmath>
#include "core_dimensional.hpp"
#include "core_dimensional_math.hpp"
#include "core_geometric.hpp"
#include "core_numeric_array.hpp"
#include "core_numeric_tuple.hpp"




//=============================================================================
namespace mhd {




//=============================================================================
using namespace dimensional;
using conserved_t                = numeric::tuple_t<unit_mass, unit_momentum, unit_momentum, unit_momentum, unit_energy>;
using conserved_density_t        = decltype(conserved_t() / unit_volume());

using flux_vector_t              = decltype(conserved_density_t() * unit_velocity());
using unit_magnetic_flux         = decltype(sqrt(unit_energy_density()) * unit_area());
using unit_magnetic_field        = decltype(sqrt(unit_energy_density()));
using unit_electric_field        = decltype(unit_magnetic_field() * unit_velocity());
using unit_vector_potential      = decltype(unit_magnetic_field() * unit_length());

using magnetic_field_vector_t    = geometric::euclidean_vector_t<unit_magnetic_field>;
using electric_field_vector_t    = geometric::euclidean_vector_t<unit_electric_field>;
using vector_potential_t         = geometric::euclidean_vector_t<unit_vector_potential>;
using velocity_vector_t          = geometric::euclidean_vector_t<unit_velocity>;

using primitive_t = numeric::tuple_t<
    unit_mass_density,
    unit_velocity,
    unit_velocity,
    unit_velocity,
    unit_energy_density,
    unit_magnetic_field,
    unit_magnetic_field,
    unit_magnetic_field>;




//=============================================================================
inline primitive_t primitive(unit_mass_density d, velocity_vector_t v, unit_energy_density p, magnetic_field_vector_t B)
{
    return numeric::tuple(
        d,
        v.component_1(),
        v.component_2(),
        v.component_3(),
        p,
        B.component_1(),
        B.component_2(),
        B.component_3());
}

inline primitive_t primitive(unit_mass_density d, unit_energy_density p)
{
    return primitive(d, {0.0, 0.0, 0.0}, p, {0.0, 0.0, 0.0});
}




//=============================================================================
inline auto mass_density    (primitive_t p) { return numeric::get<0>(p); }
inline auto velocity_1      (primitive_t p) { return numeric::get<1>(p); }
inline auto velocity_2      (primitive_t p) { return numeric::get<2>(p); }
inline auto velocity_3      (primitive_t p) { return numeric::get<3>(p); }
inline auto gas_pressure    (primitive_t p) { return numeric::get<4>(p); }
inline auto magnetic_field_1(primitive_t p) { return numeric::get<5>(p); }
inline auto magnetic_field_2(primitive_t p) { return numeric::get<6>(p); }
inline auto magnetic_field_3(primitive_t p) { return numeric::get<7>(p); }




//=============================================================================
inline auto conserved_mass      (conserved_t u) { return numeric::get<0>(u); }
inline auto conserved_momentum_1(conserved_t u) { return numeric::get<1>(u); }
inline auto conserved_momentum_2(conserved_t u) { return numeric::get<2>(u); }
inline auto conserved_momentum_3(conserved_t u) { return numeric::get<3>(u); }
inline auto conserved_energy    (conserved_t u) { return numeric::get<4>(u); }




//=============================================================================
inline auto conserved_mass_density      (conserved_density_t u) { return numeric::get<0>(u); }
inline auto conserved_momentum_density_1(conserved_density_t u) { return numeric::get<1>(u); }
inline auto conserved_momentum_density_2(conserved_density_t u) { return numeric::get<2>(u); }
inline auto conserved_momentum_density_3(conserved_density_t u) { return numeric::get<3>(u); }
inline auto conserved_energy_density    (conserved_density_t u) { return numeric::get<4>(u); }

inline auto conserved_mass_density      (flux_vector_t f) { return numeric::get<0>(f); }
inline auto conserved_momentum_density_1(flux_vector_t f) { return numeric::get<1>(f); }
inline auto conserved_momentum_density_2(flux_vector_t f) { return numeric::get<2>(f); }
inline auto conserved_momentum_density_3(flux_vector_t f) { return numeric::get<3>(f); }
inline auto conserved_energy_density    (flux_vector_t f) { return numeric::get<4>(f); }




//=============================================================================
inline auto velocity_vector(primitive_t p)
{
    return geometric::euclidean_vector(
        velocity_1(p),
        velocity_2(p),
        velocity_3(p));
}

inline auto magnetic_field_vector(primitive_t p)
{
    return geometric::euclidean_vector(
        magnetic_field_1(p),
        magnetic_field_2(p),
        magnetic_field_3(p));
}

inline auto electric_field_vector(primitive_t p)
{
    return -cross(velocity_vector(p), magnetic_field_vector(p));
}

inline auto enthalpy_density(primitive_t p, double gamma_law_index)
{
    return gas_pressure(p) * (1.0 + 1.0 / (gamma_law_index - 1.0));
}

inline auto specific_enthalpy(primitive_t p, double gamma_law_index)
{
    return enthalpy_density(p, gamma_law_index) / mass_density(p);
}

inline double specific_entropy(primitive_t p, double gamma_law_index)
{
    return std::log(gas_pressure(p).value / std::pow(mass_density(p).value, gamma_law_index));
}




//=============================================================================
inline auto sound_speed_squared(primitive_t p, double gamma_law_index)
{
    return gamma_law_index * gas_pressure(p) / mass_density(p);
}

inline auto alfven_speed_squared(primitive_t p)
{
    return length_squared(magnetic_field_vector(p)) / mass_density(p);
}

inline auto alfven_speed_squared(primitive_t p, geometric::unit_vector_t nhat)
{
    return pow<2>(dot(magnetic_field_vector(p), nhat)) / mass_density(p);
}

inline auto fast_speed_squared(primitive_t p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto cs2 = sound_speed_squared(p, gamma_law_index);
    auto cf2 = cs2 + alfven_speed_squared(p);
    auto cg2 = sqrt(cf2 * cf2 - 4.0 * cs2 * alfven_speed_squared(p, nhat));
    return 0.5 * (cf2 + cg2);
}

inline auto slow_speed_squared(primitive_t p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto cs2 = sound_speed_squared(p, gamma_law_index);
    auto cf2 = cs2 + alfven_speed_squared(p);
    auto cg2 = sqrt(cf2 * cf2 - 4.0 * cs2 * alfven_speed_squared(p, nhat));
    return 0.5 * (cf2 - cg2);
}

inline auto outer_wavespeeds(primitive_t p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto vn = dot(nhat, velocity_vector(p));
    auto cf = sqrt(fast_speed_squared(p, nhat, gamma_law_index));
    return std::pair(vn - cf, vn + cf);
}




//=============================================================================
inline auto momentum_density_vector(conserved_density_t u)
{
    return geometric::euclidean_vector(
        conserved_momentum_density_1(u),
        conserved_momentum_density_2(u),
        conserved_momentum_density_3(u));
}

inline auto momentum_density_vector(flux_vector_t f)
{
    return geometric::euclidean_vector(
        conserved_momentum_density_1(f),
        conserved_momentum_density_2(f),
        conserved_momentum_density_3(f));
}

inline auto kinetic_energy_density(conserved_density_t u)
{
    return 0.5 * length_squared(momentum_density_vector(u)) / conserved_mass_density(u);
}

inline auto magnetic_energy_density(magnetic_field_vector_t b)
{
    return 0.5 * length_squared(b);
}

inline auto internal_energy_density(conserved_density_t u, magnetic_field_vector_t b)
{
    return conserved_energy_density(u) - kinetic_energy_density(u) - magnetic_energy_density(b);
}




//=============================================================================
inline primitive_t recover_primitive(conserved_density_t u, magnetic_field_vector_t b, double gamma_law_index)
{
    auto d = conserved_mass_density(u);
    auto v = momentum_density_vector(u) / d;
    auto p = internal_energy_density(u, b) * (gamma_law_index - 1);
    return primitive(d, v, p, b);
}

inline conserved_density_t conserved_density(primitive_t p, double gamma_law_index)
{
    auto d = mass_density(p);
    auto v = velocity_vector(p);
    auto B = magnetic_field_vector(p);
    auto E = 0.5 * d * length_squared(v) + gas_pressure(p) / (gamma_law_index - 1.0) + magnetic_energy_density(B);

    return numeric::tuple(
        d,
        d * v.component_1(),
        d * v.component_2(),
        d * v.component_3(),
        E);
}




//=============================================================================
inline flux_vector_t flux(primitive_t p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto v = velocity_vector(p);
    auto B = magnetic_field_vector(p);
    auto Bn = dot(B, nhat);
    auto vn = dot(v, nhat);
    auto pt = gas_pressure(p) + 0.5 * magnetic_energy_density(B);

    auto pressure_term = numeric::tuple(
        unit_mass_density() * unit_velocity(),
        pt * nhat.component_1(),
        pt * nhat.component_2(),
        pt * nhat.component_3(),
        pt * vn);

    auto maxwell_term = numeric::tuple(
        unit_mass_density() * unit_velocity(),
        Bn * B.component_1(),
        Bn * B.component_2(),
        Bn * B.component_3(),
        Bn * dot(B, v));

    return vn * conserved_density(p, gamma_law_index) + pressure_term + maxwell_term;
}

inline electric_field_vector_t induction(primitive_t p, geometric::unit_vector_t nhat)
{
    return cross(nhat, electric_field_vector(p));
}




//=============================================================================
inline auto riemann_hlle(
    primitive_t pl,
    primitive_t pr,
    unit_magnetic_field b_longitudinal,
    geometric::unit_vector_t nhat,
    double gamma_law_index)
{
    auto ul = conserved_density(pl, gamma_law_index);
    auto ur = conserved_density(pr, gamma_law_index);
    auto fl = flux(pl, nhat, gamma_law_index);
    auto fr = flux(pr, nhat, gamma_law_index);

    auto bl = magnetic_field_vector(pl);
    auto br = magnetic_field_vector(pr);
    auto il = induction(pl, nhat);
    auto ir = induction(pr, nhat);

    auto [alm, alp] = outer_wavespeeds(pl, nhat, gamma_law_index);
    auto [arm, arp] = outer_wavespeeds(pr, nhat, gamma_law_index);
    auto ap = std::max(unit_velocity(0), std::max(alp, arp));
    auto am = std::min(unit_velocity(0), std::min(alm, arm));

    auto hll_flux      = (ap * fl - am * fr - (ul - ur) * ap * am) / (ap - am);
    auto hll_induction = (ap * il - am * ir - (bl - br) * ap * am) / (ap - am);

    // -n x F(B) = -n x (n x E) = -n (n.E) + E (n.n) = -E_parallel + E = E_trans
    auto hll_electric_field = -cross(nhat, hll_induction);

    return std::pair(hll_flux, hll_electric_field);
}

} // namespace mhd




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_mhd()
{
    using namespace mhd;

    auto require_cons_to_prim = [] (mhd::primitive_t p0)
    {
        auto u0 = conserved_density(p0, 5. / 3);
        auto B0 = magnetic_field_vector(p0);
        auto p1 = recover_primitive(u0, B0, 5. / 3);
        require(all(map(p1 - p0, [] (auto x) { return std::abs(x.value) < 1e-10; })));
    };

    require_cons_to_prim(primitive(1.0, 1.0));
    require_cons_to_prim(primitive(1.5, 1e3));
    require_cons_to_prim(primitive(1.0, {0.0, 0.0, 0.0}, 1.0, {0.0, 0.0, 0.0}));
    require_cons_to_prim(primitive(5.7, {1.0, 2.0, 3.0}, 5.0, {0.0, 0.0, 0.0}));
    require_cons_to_prim(primitive(5.7, {1.0, 2.0, 3.0}, 5.0, {3.0, 2.0, 1.0}));

    auto n1 = geometric::unit_vector_on(1);
    auto n2 = geometric::unit_vector_on(2);
    auto n3 = geometric::unit_vector_on(3);
    auto p = mhd::primitive(1.0, {1.0, 0.0, 0.0}, 1.0, {0.0, 1.0, 0.0});
    require(-cross(n1, mhd::induction(p, n1)) == electric_field_vector(p));
    require(-cross(n2, mhd::induction(p, n2)) == electric_field_vector(p));
    require(-cross(n3, mhd::induction(p, n3)) == electric_field_vector(p) * 0.0);
}

#endif // DO_UNIT_TESTS
