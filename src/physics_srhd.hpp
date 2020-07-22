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
namespace srhd {


//=============================================================================
using namespace dimensional;
using primitive_t = numeric::tuple_t<unit_mass_density, unit_scalar, unit_scalar, unit_scalar, unit_energy_density>;
using conserved_t = numeric::tuple_t<unit_mass, unit_momentum, unit_momentum, unit_momentum, unit_energy>;
using conserved_density_t = decltype(conserved_t() / unit_volume());
using flux_vector_t       = decltype(conserved_density_t() * unit_velocity());
const static auto light_speed = dimensional::unit_velocity(1.0);




//=============================================================================
template<typename ExceptionType>
#ifdef SRHD_NO_EXCEPTIONS
void throw_or_abort(std::string message) noexcept
{
    std::printf("FATAL: %s\n", message.data());
    exit(1);
}
#else
#include <stdexcept>
void throw_or_abort(std::string message)
{
    throw ExceptionType(message);
}
#endif




//=============================================================================
inline primitive_t primitive(unit_mass_density d, unit_scalar u, unit_scalar v, unit_scalar w, unit_energy_density p)
{
    return {{d, u, v, w, p}};
}

inline primitive_t primitive(unit_mass_density d, unit_energy_density p)
{
    return primitive(d, 0.0, 0.0, 0.0, p);
}




//=============================================================================
inline auto mass_density(primitive_t p) { return numeric::get<0>(p); }
inline auto gamma_beta_1(primitive_t p) { return numeric::get<1>(p); }
inline auto gamma_beta_2(primitive_t p) { return numeric::get<2>(p); }
inline auto gamma_beta_3(primitive_t p) { return numeric::get<3>(p); }
inline auto gas_pressure(primitive_t p) { return numeric::get<4>(p); }


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
inline auto gamma_beta_squared(primitive_t p)
{
    auto u1 = gamma_beta_1(p);
    auto u2 = gamma_beta_2(p);
    auto u3 = gamma_beta_3(p);
    return u1 * u1 + u2 * u2 + u3 * u3;
}

inline auto lorentz_factor_squared(primitive_t p)
{
    return 1.0 + gamma_beta_squared(p);
}

inline auto lorentz_factor(primitive_t p)
{
    return std::sqrt(lorentz_factor_squared(p));
}




//=============================================================================
inline auto velocity_1(primitive_t p) { return light_speed * gamma_beta_1(p) / lorentz_factor(p); }
inline auto velocity_2(primitive_t p) { return light_speed * gamma_beta_2(p) / lorentz_factor(p); }
inline auto velocity_3(primitive_t p) { return light_speed * gamma_beta_3(p) / lorentz_factor(p); }
inline auto velocity_vector(primitive_t p)
{
    auto gamma = lorentz_factor(p);

    return geometric::euclidean_vector(
        light_speed * gamma_beta_1(p) / gamma,
        light_speed * gamma_beta_2(p) / gamma,
        light_speed * gamma_beta_3(p) / gamma);
}




//=============================================================================
inline auto enthalpy_density(primitive_t p, double gamma_law_index)
{
    auto c2 = light_speed * light_speed;
    return mass_density(p) * c2 + gas_pressure(p) * (1.0 + 1.0 / (gamma_law_index - 1.0));
}

inline auto specific_enthalpy(primitive_t p, double gamma_law_index)
{
    return enthalpy_density(p, gamma_law_index) / mass_density(p);
}

inline double specific_entropy(primitive_t p, double gamma_law_index)
{
    return std::log(gas_pressure(p).value / std::pow(mass_density(p).value, gamma_law_index));
}

inline auto sound_speed_squared(primitive_t p, double gamma_law_index)
{
    auto c2 = light_speed * light_speed;
    return gamma_law_index * c2 * gas_pressure(p) / enthalpy_density(p, gamma_law_index);
}

inline auto outer_wavespeeds(primitive_t p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto c2 = light_speed * light_speed;
    auto a2 = sound_speed_squared(p, gamma_law_index) / c2;
    auto uu = gamma_beta_squared(p);
    auto vn = dot(nhat, velocity_vector(p));
    auto bn = vn / light_speed;
    auto bb = uu / (1 + uu);
    auto b2 = bn * bn;
    auto k0 = sqrt(a2 * (1 - bb) * (1 - bb * a2 - b2 * (1 - a2))) * light_speed;

    return std::pair((vn * (1 - a2) - k0) / (1 - bb * a2),
                     (vn * (1 - a2) + k0) / (1 - bb * a2));
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

inline auto internal_energy_density(conserved_density_t u)
{
    return conserved_energy_density(u) - kinetic_energy_density(u);
}




//=============================================================================
inline primitive_t recover_primitive(conserved_density_t u, double gamma_law_index, double temperature_floor=0.0)
{
    auto newton_iter_max = 50;
    auto error_tolerance = 1e-12;
    auto c2        = light_speed * light_speed;
    auto gm        = gamma_law_index;
    auto D         = conserved_mass_density(u);
    auto tau       = conserved_energy_density(u);
    auto SS        = length_squared(momentum_density_vector(u));
    auto found     = false;
    auto iteration = 0;
    auto W0        = 1.0;
    auto p         = dimensional::unit_energy_density(0.0);

    while (iteration <= newton_iter_max)
    {
        auto Et = tau + p + D * c2;
        auto b2 = std::min(SS * c2 / Et / Et, unit_scalar(1.0 - 1e-10));
        auto W2 = 1.0 / (1.0 - b2);
        auto W  = std::sqrt(W2);
        auto e  = (tau + D * c2 * (1.0 - W) + p * (1.0 - W2)) / (D * W);
        auto d  = D / W;
        auto h  = c2 + e + p / d;
        auto a2 = gm * p / (d * h);
        auto f  = d * e * (gm - 1.0) - p;
        auto g  = b2 * a2 - 1.0;

        p -= f / g;

        if (std::fabs(f.value) < error_tolerance || (
            std::fabs(f.value) < error_tolerance * 1e4 && iteration == newton_iter_max))
        {
            W0 = W;
            found = true;
            break;
        }
        ++iteration;
    }

    if (temperature_floor > 0.0)
        p = std::max(p, temperature_floor * D * c2);

    if (! found)
        throw_or_abort<std::invalid_argument>("srhd::recover_primitive (root finder not converging)");

    if (temperature_floor <= 0.0 && p <= unit_energy_density(0.0))
        throw_or_abort<std::invalid_argument>("srhd::recover_primitive (negative pressure)");

    return numeric::tuple(
        D / W0,
        W0 * conserved_momentum_density_1(u) / (tau + D * c2 + p) * light_speed,
        W0 * conserved_momentum_density_2(u) / (tau + D * c2 + p) * light_speed,
        W0 * conserved_momentum_density_3(u) / (tau + D * c2 + p) * light_speed,
        p);
}




//=============================================================================
inline conserved_density_t conserved_density(primitive_t p, double gamma_law_index)
{
    auto c = light_speed;
    auto W = lorentz_factor(p);
    auto D = mass_density(p) * W;
    auto h = specific_enthalpy(p, gamma_law_index);

    return numeric::tuple(
        D,
        D * h * gamma_beta_1(p) / c,
        D * h * gamma_beta_2(p) / c,
        D * h * gamma_beta_3(p) / c,
        D * (h * W - c * c) - gas_pressure(p));
}




//=============================================================================
inline flux_vector_t flux(primitive_t p, geometric::unit_vector_t nhat, double gamma_law_index)
{
    auto pg = gas_pressure(p);
    auto vn = dot(nhat, velocity_vector(p));
    auto md = unit_mass_density(0.0) * unit_velocity(0.0);

    auto pressure_term = numeric::tuple(
        md,
        pg * nhat.component_1(),
        pg * nhat.component_2(),
        pg * nhat.component_3(),
        pg * vn);

    return vn * conserved_density(p, gamma_law_index) + pressure_term;
}




//=============================================================================
inline auto spherical_geometry_source_terms(primitive_t p, unit_length spherical_radius, unit_scalar polar_angle_theta, double gamma_law_index)
{
    auto cotq = std::tan(M_PI_2 - polar_angle_theta);
    auto ur = gamma_beta_1(p);
    auto uq = gamma_beta_2(p);
    auto up = gamma_beta_3(p);
    auto pg = gas_pressure(p);
    auto H0 = enthalpy_density(p, gamma_law_index);
    auto Sd = unit_mass_density(0.0) / unit_time(1.0);
    auto Sr = (2.0  * pg + H0 * (uq * uq        + up * up)) / spherical_radius;
    auto Sq = (cotq * pg + H0 * (up * up * cotq - ur * uq)) / spherical_radius;
    auto Sp =        -up * H0 * (ur + uq * cotq) / spherical_radius;
    auto Se = unit_energy_density(0.0) / unit_time(1.0);

    return numeric::tuple(Sd, Sr, Sq, Sp, Se);
}




//=============================================================================
struct riemann_solver_mode_hlle_fluxes_t {};
struct riemann_solver_mode_hllc_fluxes_t {};
struct riemann_solver_mode_hlle_fluxes_moving_face_t { unit_velocity face_speed; };
struct riemann_solver_mode_hllc_fluxes_moving_face_t { unit_velocity face_speed; };
struct riemann_solver_mode_hllc_fluxes_across_contact_t {};
struct riemann_solver_mode_contact_speed_t {};




/**
 * @brief      HLLC Riemann solver based on Mignone & Bodo (2005)
 *
 * @param[in]  pl                 Primitive state left of the interface
 * @param[in]  pr                 Primitive state right of the interface
 * @param[in]  nhat               The unit vector normal to the interface
 * @param[in]  gamma_law_index    The adiabatic gamma law index
 * @param[in]  mode               Instance of one of the structs above,
 *                                determining what will be returned
 *
 * @tparam     RiemannSolverMode  The type of the calculation mode
 *
 * @return     Either the contact speed, or fluxes, depending on the mode
 *
 * @note       Equation numbers given are from Mignone & Bodo (2005).
 */
template<typename RiemannSolverMode>
auto riemann_solver(primitive_t pl, primitive_t pr, geometric::unit_vector_t nhat, double gamma_law_index, RiemannSolverMode mode)
{
    auto ul = conserved_density(pl, gamma_law_index);
    auto ur = conserved_density(pr, gamma_law_index);
    auto fl = flux(pl, nhat, gamma_law_index);
    auto fr = flux(pr, nhat, gamma_law_index);

    auto [alm, alp] = outer_wavespeeds(pl, nhat, gamma_law_index);
    auto [arm, arp] = outer_wavespeeds(pr, nhat, gamma_law_index);
    auto ar = std::max(alp, arp);
    auto al = std::min(alm, arm);

    if constexpr (std::is_same_v<RiemannSolverMode, riemann_solver_mode_hlle_fluxes_t>)
    {
        if (unit_velocity(0.0) < al) return fl;
        if (unit_velocity(0.0) > ar) return fr;
        return (ar * fl - al * fr - (ul - ur) * ar * al) / (ar - al);
    } else {

    // Equations (9) and (11)
    auto u_hll = (ar * ur - al * ul + (fl - fr))           / (ar - al);
    auto f_hll = (ar * fl - al * fr - (ul - ur) * ar * al) / (ar - al);

    if constexpr (std::is_same_v<RiemannSolverMode, riemann_solver_mode_hlle_fluxes_moving_face_t>)
    {
        if (mode.face_speed < al) return fl - mode.face_speed * ul;
        if (mode.face_speed > ar) return fr - mode.face_speed * ur;
        return f_hll - mode.face_speed * u_hll;
    } else {

    auto c2 = light_speed * light_speed;

    // Mignone defines total energy to include rest mass
    auto ue_hll = conserved_energy_density(u_hll) + conserved_mass_density(u_hll) * c2;
    auto fe_hll = conserved_energy_density(f_hll) + conserved_mass_density(f_hll) * c2;
    auto um_hll = dot(momentum_density_vector(u_hll), nhat);
    auto fm_hll = dot(momentum_density_vector(f_hll), nhat);

    // Equation (18)
    auto a = +fe_hll;
    auto b = -fm_hll - ue_hll;
    auto c = +um_hll;
    auto d = b * b - 4.0 * a * c;
    auto as = d.value < 0.0 ? unit_velocity(0.0) : (std::abs(a.value) > 1e-8 ? (-b - sqrt(d)) / (2.0 * a) : -c / b) * c2;

    if (std::isnan(as.value) || as.value < -1e3 || as.value > 1e3)
    {
        std::printf("pl = %e %e %e\n", numeric::get<0>(pl).value, numeric::get<1>(pl).value, numeric::get<4>(pl).value);
        std::printf("pr = %e %e %e\n", numeric::get<0>(pr).value, numeric::get<1>(pr).value, numeric::get<4>(pr).value);

        std::printf("a = %lf\n", a.value);
        std::printf("b = %lf\n", b.value);
        std::printf("b^2 - 4 a c = %lf\n", d.value);
        throw std::invalid_argument("srhd::riemann_solver (nan contact speed)");
    }
    if constexpr (std::is_same_v<RiemannSolverMode, riemann_solver_mode_contact_speed_t>)
    {
        return as;
    } else {

    auto ps = -fe_hll * as / c2 + fm_hll;

    auto star_state_flux = [nhat, c2, as, ps] (conserved_density_t u, flux_vector_t f, primitive_t p, unit_velocity a, unit_velocity face_speed)
    {
        // Equations (16)
        auto S = momentum_density_vector(u);
        auto v = dot(velocity_vector(p), nhat);
        auto e = conserved_energy_density(u) + conserved_mass_density(u) * c2;
        auto m = dot(S, nhat);
        auto n = S - m * nhat; // transverse momentum
        auto es = (e * (a - v) + ps * as - gas_pressure(p) * v) / (a - as);
        auto ms = (m * (a - v) + ps - gas_pressure(p))          / (a - as);
        auto ns =  n * (a - v)                                  / (a - as);
        auto Ds = conserved_mass_density(u) * (a - v)           / (a - as);
        auto Ss = ms * nhat + ns;

        auto us = numeric::tuple(Ds, Ss.component_1(), Ss.component_2(), Ss.component_3(), es - Ds * c2);
        auto fs = f + a * (us - u);

        return fs - face_speed * us;
    };

    if constexpr (std::is_same_v<RiemannSolverMode, riemann_solver_mode_hllc_fluxes_across_contact_t>)
    {
        return std::pair(star_state_flux(ul, fl, pl, al, as), as);
    }

    if constexpr (std::is_same_v<RiemannSolverMode, riemann_solver_mode_hllc_fluxes_t>)
    {
        if (unit_velocity(0.0) <= al) return fl;
        if (unit_velocity(0.0) >= ar) return fr;
        if (unit_velocity(0.0) <= as) return star_state_flux(ul, fl, pl, al, 0.0);
        if (unit_velocity(0.0) >= as) return star_state_flux(ur, fr, pr, ar, 0.0);

        return flux_vector_t{};        
    }

    if constexpr (std::is_same_v<RiemannSolverMode, riemann_solver_mode_hllc_fluxes_moving_face_t>)
    {
        if (mode.face_speed <= al) return fl - mode.face_speed * ul;
        if (mode.face_speed >= ar) return fr - mode.face_speed * ur;
        if (mode.face_speed <= as) return star_state_flux(ul, fl, pl, al, mode.face_speed);
        if (mode.face_speed >= as) return star_state_flux(ur, fr, pr, ar, mode.face_speed);

        return flux_vector_t{};
    }
    }
    }
    }
}

inline flux_vector_t riemann_hlle(primitive_t pl, primitive_t pr, geometric::unit_vector_t nhat, double gamma_law_index)
{
    return riemann_solver(pl, pr, nhat, gamma_law_index, riemann_solver_mode_hlle_fluxes_t{});
}

inline flux_vector_t riemann_hlle(primitive_t pl, primitive_t pr, geometric::unit_vector_t nhat, unit_velocity face_speed, double gamma_law_index)
{
    return riemann_solver(pl, pr, nhat, gamma_law_index, riemann_solver_mode_hlle_fluxes_moving_face_t{face_speed});
}

inline flux_vector_t riemann_hllc(primitive_t pl, primitive_t pr, geometric::unit_vector_t nhat, double gamma_law_index)
{
    return riemann_solver(pl, pr, nhat, gamma_law_index, riemann_solver_mode_hllc_fluxes_t{});
}

inline flux_vector_t riemann_hllc(primitive_t pl, primitive_t pr, geometric::unit_vector_t nhat, unit_velocity face_speed, double gamma_law_index)
{
    return riemann_solver(pl, pr, nhat, gamma_law_index, riemann_solver_mode_hllc_fluxes_moving_face_t{face_speed});
}

} // namespace srhd




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_srhd()
{
    using namespace srhd;

    auto require_cons_to_prim = [] (srhd::primitive_t p0)
    {
        auto p1 = recover_primitive(conserved_density(p0, 4. / 3), 4. / 3);
        require(all(map(p1 - p0, [] (auto x) { return std::abs(x.value) < 1e-10; })));
    };

    require_cons_to_prim(primitive(1.0, 1.0));
    require_cons_to_prim(primitive(1.5, 1e3));
    require_cons_to_prim(primitive(1.0, 0.0, 0.0, 0.0, 1.0));
    require_cons_to_prim(primitive(5.7, 1.0, 2.0, 3.0, 5.0));

    auto require_contact_speed = [] (unsigned dir)
    {
        auto xh = geometric::unit_vector_on(dir);
        auto pl = primitive(1.0, 0.5 * (dir == 1), 0.6 * (dir == 2), 0.7 * (dir == 3), 10.0);
        auto pr = primitive(4.0, 0.5 * (dir == 1), 0.6 * (dir == 2), 0.7 * (dir == 3), 10.0);
        auto contact_speed = riemann_solver(pl, pr, xh, 4. / 3, riemann_solver_mode_contact_speed_t());
        auto fc = riemann_solver(pl, pr, xh, 4. / 3, riemann_solver_mode_hllc_fluxes_moving_face_t{contact_speed});

        require(abs(contact_speed - velocity_vector(pl).component(dir)).value < 1e-12);
        require(abs(conserved_mass_density(fc)).value < 1e-12);
    };
    require_contact_speed(1);
    require_contact_speed(2);
    require_contact_speed(3);
}

#endif // DO_UNIT_TESTS
