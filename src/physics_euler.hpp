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
namespace euler {


//=============================================================================
using namespace dimensional;
using primitive_t = numeric::tuple_t<unit_mass_density, unit_velocity, unit_velocity, unit_velocity, unit_energy_density>;
using conserved_t = numeric::tuple_t<unit_mass, unit_momentum, unit_momentum, unit_momentum, unit_energy>;
using conserved_density_t = decltype(conserved_t{} / unit_volume{});
using flux_vector_t       = decltype(conserved_density_t{} * unit_velocity{});




//=============================================================================
inline primitive_t primitive(unit_mass_density d, unit_velocity u, unit_velocity v, unit_velocity w, unit_energy_density p)
{
    return {{d, u, v, w, p}};
}

inline primitive_t primitive(unit_mass_density d, unit_energy_density p)
{
    return primitive(d, 0.0, 0.0, 0.0, p);
}




//=============================================================================
inline const auto& mass_density(const primitive_t& p) { return numeric::get<0>(p); }
inline const auto& velocity_1  (const primitive_t& p) { return numeric::get<1>(p); }
inline const auto& velocity_2  (const primitive_t& p) { return numeric::get<2>(p); }
inline const auto& velocity_3  (const primitive_t& p) { return numeric::get<3>(p); }
inline const auto& gas_pressure(const primitive_t& p) { return numeric::get<4>(p); }


//=============================================================================
inline const auto& conserved_mass      (const conserved_t& u) { return numeric::get<0>(u); }
inline const auto& conserved_momentum_1(const conserved_t& u) { return numeric::get<1>(u); }
inline const auto& conserved_momentum_2(const conserved_t& u) { return numeric::get<2>(u); }
inline const auto& conserved_momentum_3(const conserved_t& u) { return numeric::get<3>(u); }
inline const auto& conserved_energy    (const conserved_t& u) { return numeric::get<4>(u); }


//=============================================================================
inline const auto& conserved_mass_density      (const conserved_density_t& u) { return numeric::get<0>(u); }
inline const auto& conserved_momentum_density_1(const conserved_density_t& u) { return numeric::get<1>(u); }
inline const auto& conserved_momentum_density_2(const conserved_density_t& u) { return numeric::get<2>(u); }
inline const auto& conserved_momentum_density_3(const conserved_density_t& u) { return numeric::get<3>(u); }
inline const auto& conserved_energy_density    (const conserved_density_t& u) { return numeric::get<4>(u); }




//=============================================================================
inline auto velocity_vector(const primitive_t& p)
{
    return geometric::euclidean_vector(
        velocity_1(p),
        velocity_2(p),
        velocity_3(p));
}

inline auto sound_speed_squared(const primitive_t& p, double gamma_law_index)
{
    return gamma_law_index * gas_pressure(p) / mass_density(p);
}

inline auto enthalpy_density(const primitive_t& p, double gamma_law_index)
{
    return gas_pressure(p) * (1.0 + 1.0 / (gamma_law_index - 1.0));
}

inline auto specific_enthalpy(const primitive_t& p, double gamma_law_index)
{
    return enthalpy_density(p, gamma_law_index) / mass_density(p);
}

inline double specific_entropy(const primitive_t& p, double gamma_law_index)
{
    return std::log(gas_pressure(p).value / std::pow(mass_density(p).value, gamma_law_index));
}

inline auto outer_wavespeeds(const primitive_t& p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto vn = dot(nhat, velocity_vector(p));
    auto cs = sqrt(sound_speed_squared(p, gamma_law_index));
    return std::pair(vn - cs, vn + cs);
}




//=============================================================================
inline auto momentum_vector(const conserved_t& U)
{
    return geometric::euclidean_vector(
        conserved_momentum_1(U),
        conserved_momentum_2(U),
        conserved_momentum_3(U));
}

inline auto momentum_density_vector(const conserved_density_t& u)
{
    return geometric::euclidean_vector(
        conserved_momentum_density_1(u),
        conserved_momentum_density_2(u),
        conserved_momentum_density_3(u));
}

inline auto kinetic_energy_density(const conserved_density_t& u)
{
    return 0.5 * length_squared(momentum_density_vector(u)) / conserved_mass_density(u);
}

inline auto internal_energy_density(const conserved_density_t& u)
{
    return conserved_energy_density(u) - kinetic_energy_density(u);
}




//=============================================================================
inline primitive_t recover_primitive(conserved_density_t u, double gamma_law_index)
{
    return numeric::tuple(
        conserved_mass_density(u),
        conserved_momentum_density_1(u) / conserved_mass_density(u),
        conserved_momentum_density_2(u) / conserved_mass_density(u),
        conserved_momentum_density_3(u) / conserved_mass_density(u),
        internal_energy_density(u) * (gamma_law_index - 1.0));
}




//=============================================================================
inline conserved_density_t conserved_density(primitive_t p, double gamma_law_index)
{
    return numeric::tuple(
        mass_density(p),
        mass_density(p) * velocity_1(p),
        mass_density(p) * velocity_2(p),
        mass_density(p) * velocity_3(p),
        mass_density(p) * length_squared(velocity_vector(p)) * 0.5 + gas_pressure(p) / (gamma_law_index - 1.0));
}




//=============================================================================
inline flux_vector_t flux(primitive_t p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto pg = gas_pressure(p);
    auto vn = dot(nhat, velocity_vector(p));
    auto md = unit_mass_density{} * unit_velocity{};

    auto pressure_term = numeric::tuple(
        md,
        pg * nhat.component_1(),
        pg * nhat.component_2(),
        pg * nhat.component_3(),
        pg * vn);

    return vn * conserved_density(p, gamma_law_index) + pressure_term;
}




//=============================================================================
inline flux_vector_t riemann_hlle( const primitive_t& pl, const primitive_t& pr, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto ul = conserved_density(pl, gamma_law_index);
    auto ur = conserved_density(pr, gamma_law_index);
    auto fl = flux(pl, nhat, gamma_law_index);
    auto fr = flux(pr, nhat, gamma_law_index);

    auto [alm, alp] = outer_wavespeeds(pl, nhat, gamma_law_index);
    auto [arm, arp] = outer_wavespeeds(pr, nhat, gamma_law_index);
    auto ap = std::max(unit_velocity{}, std::max(alp, arp));
    auto am = std::min(unit_velocity{}, std::min(alm, arm));

    return (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am);
}

} // namespace euler
