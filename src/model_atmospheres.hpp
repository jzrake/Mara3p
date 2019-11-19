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
#include "core_numeric_tuple.hpp"




//=============================================================================
namespace mara {

using namespace dimensional;




//=============================================================================
namespace detail
{

template<typename Function, typename T, typename U>
T solve_secant(Function f, T starting_guess1, T starting_guess2, U tolerance)
{
    T x1 = starting_guess1;
    T x2 = starting_guess2;
    auto y1 = f(x1);
    auto y2 = f(x2);

    while (abs(y2) > tolerance)
    {
        auto x_next = x2 - y2 * (x2 - x1) / (y2 - y1);
        auto y_next = f(x_next);

        x1 = x2;
        y1 = y2;
        x2 = x_next;
        y2 = y_next;
    }
    return x2;
}

} // namespace detail




//=============================================================================
struct cold_power_law_medium_t
{


    //=========================================================================
    cold_power_law_medium_t with_inner_radius(unit_length new_inner_radius) const
    {
        return {new_inner_radius, density_at_base, density_index};
    }

    cold_power_law_medium_t with_density_at_base(unit_mass_density new_density_at_base) const
    {
        return {inner_radius, new_density_at_base, density_index};
    }

    cold_power_law_medium_t with_density_index(double new_density_index) const
    {
        return {inner_radius, density_at_base, new_density_index};
    }


    //=========================================================================
    unit_mass_density density_at(unit_length r, unit_time t) const
    {
        return density_at_base * std::pow(r / inner_radius, -density_index);
    }

    auto primitive_srhd(unit_length r, unit_time t, unit_scalar entropy=-6.0) const
    {
        auto d = density_at(r, t);
        auto p = unit_energy_density(std::pow(d.value, 4. / 3) * std::exp(entropy));
        return numeric::tuple(d, unit_scalar(0.0), unit_scalar(0.0), unit_scalar(0.0), p);
    }


    //=========================================================================
    unit_length       inner_radius    = 1.0;
    unit_mass_density density_at_base = 1.0;
    unit_velocity     light_speed     = 1.0;
    double            density_index   = 2.0;
};




//=============================================================================
struct cold_wind_model_t
{


    //=========================================================================
    cold_wind_model_t with_kinetic_luminosity(unit_power new_luminosity) const
    {
        return {new_luminosity, time_scale, gamma_beta0, solid_angle};
    }

    cold_wind_model_t with_time_scale(unit_time new_time_scale) const
    {
        return {luminosity, new_time_scale, gamma_beta0, solid_angle};
    }

    cold_wind_model_t with_gamma_beta(unit_scalar new_gamma_beta) const
    {
        return {luminosity, time_scale, new_gamma_beta, solid_angle};
    }

    cold_wind_model_t with_solid_angle(unit_scalar new_solid_angle) const
    {
        return {luminosity, time_scale, gamma_beta0, new_solid_angle};
    }


    //=========================================================================
    unit_energy total_energy() const
    {
        auto gu0 = lorentz_factor(0) * gamma_beta(0);
        return 2. / 3 * luminosity / gu0 * time_scale * (pow<3, 2>(unit_scalar(1.0) + gamma_beta0) - 1.0);
    }

    unit_scalar lorentz_factor(unit_time t) const
    {
        return std::sqrt(1.0 + pow<2>(gamma_beta(t)));
    }

    unit_scalar gamma_beta(unit_time t) const
    {
        return gamma_beta0 * std::exp(-(t - unit_time(1.0)) / time_scale);
    }

    unit_power kinetic_luminosity(unit_time t) const
    {
        auto gu0 = lorentz_factor(0) * gamma_beta(0);
        auto gu1 = lorentz_factor(t) * gamma_beta(t);
        return luminosity * gu1 / gu0;
    }

    unit_mass_rate mass_loss_rate(unit_time t) const
    {
        return kinetic_luminosity(t) / lorentz_factor(t) / pow<2>(light_speed);
    }


    //=========================================================================
    unit_mass_density density_at(unit_length r, unit_time t) const
    {
        return mass_loss_rate(t) / (solid_angle * r * r * gamma_beta(t) * light_speed);
    }

    auto primitive_srhd(unit_length r, unit_time t, unit_scalar entropy=-6.0) const
    {
        auto d = density_at(r, t);
        auto u = gamma_beta(t);
        auto p = unit_energy_density(std::pow(d.value, 4. / 3) * std::exp(entropy));
        return numeric::tuple(d, u, unit_scalar(0.0), unit_scalar(0.0), p);
    }


    //=========================================================================
    unit_power    luminosity    = 1.0; // = gamma rho u c^3 solid_angle r^2;
    unit_time     time_scale    = 1.0;
    unit_scalar   gamma_beta0   = 1.0;
    unit_scalar   solid_angle   = 4 * M_PI;
    unit_velocity light_speed   = 1.0;
};




//=============================================================================
struct time_varying_cold_wind_t
{


    //=========================================================================
    using mass_loss_rate_function_t = std::function<unit_mass_rate(unit_time)>;
    using gamma_beta_function_t     = std::function<unit_scalar(unit_time)>;


    //=========================================================================
    auto with_solid_angle(unit_scalar new_solid_angle) const
    -> time_varying_cold_wind_t
    {
        return {new_solid_angle, light_speed, mass_loss_rate, gamma_beta};
    }

    auto with_mass_loss_rate(mass_loss_rate_function_t new_mass_loss_rate) const
    -> time_varying_cold_wind_t
    {
        return {solid_angle, light_speed, new_mass_loss_rate, gamma_beta};
    }

    auto with_gamma_beta(gamma_beta_function_t new_gamma_beta) const
    -> time_varying_cold_wind_t
    {
        return {solid_angle, light_speed, mass_loss_rate, new_gamma_beta};
    }


    //=========================================================================
    unit_scalar lorentz_factor(unit_time t) const
    {
        return std::sqrt(1.0 + pow<2>(gamma_beta(t)));
    }

    unit_power kinetic_luminosity(unit_time t) const
    {
        return lorentz_factor(t) * mass_loss_rate(t) * pow<2>(light_speed);
    }


    //=========================================================================
    unit_mass_density density_at(unit_length r, unit_time t) const
    {
        return mass_loss_rate(t) / (solid_angle * r * r * gamma_beta(t) * light_speed);
    }

    auto primitive_srhd(unit_length r, unit_time t, unit_scalar log10_temperature=-6.0) const
    {
        auto d = density_at(r, t);
        auto u = gamma_beta(t);
        auto p = d * pow<2>(light_speed) * std::pow(10.0, log10_temperature);
        return numeric::tuple(d, u, unit_scalar(0.0), unit_scalar(0.0), p);
    }


    //=========================================================================
    unit_scalar   solid_angle = 4 * M_PI;
    unit_velocity light_speed = 1.0;
    mass_loss_rate_function_t mass_loss_rate;
    gamma_beta_function_t     gamma_beta;
};




/**
 * @brief      A model describing an expanding gas cloud, with power-law density
 *             profile, surrounded by a relativistic envelop. The envelop is an
 *             ensemble of balistic shells ejected simultaneously from the
 *             origin at t=0. The shell with mass dm contains gas with
 *             four-velocity between u and u + du. Thus the envelop profile is
 *             given by the distribution du / dm or equivalently by u(m) = u1 *
 *             (m / m1)^(-psi). Here u1 is the four-velocity of the fastest
 *             shell, and m1 is the mass outside the fastest shell. Note that
 *             u(m) is a monotonically decreasing function of mass, so the
 *             slowest-moving shell is determined by the total envelop mass,
 *             nominally ~ 0.5% of a solar mass. The cloud is a cold wind whose
 *             density and velocity match that of the slowest moving shell at r
 *             = cloud_outer_boundary(t) = v_cloud * t. Thus u_cloud =
 *             u(m_envelop) and the cloud density decreases as r^-2.
 */
struct cloud_with_envelop_model_t
{


    //=========================================================================
    cloud_with_envelop_model_t with_inner_radius(unit_length r0) const
    {
        auto result = *this;
        result.inner_radius = r0;
        return result;
    }

    cloud_with_envelop_model_t with_envelop_mass(unit_mass m) const
    {
        auto result = *this;
        result.envelop_mass = m;
        return result;
    }

    cloud_with_envelop_model_t with_fastest_shell_mass(unit_mass m) const
    {
        auto result = *this;
        result.fastest_shell_mass = m;
        return result;
    }

    cloud_with_envelop_model_t with_fastest_gamma_beta(unit_scalar u) const
    {
        auto result = *this;
        result.fastest_gamma_beta = u;
        return result;
    }

    cloud_with_envelop_model_t with_solid_angle(unit_scalar omega) const
    {
        auto result = *this;
        result.solid_angle = omega;
        return result;        
    }

    cloud_with_envelop_model_t with_psi(double psi) const
    {
        auto result = *this;
        result.psi = psi;
        return result;        
    }


    //=========================================================================
    unit_scalar gamma_beta(unit_mass m) const
    {
        return fastest_gamma_beta * std::pow(m / fastest_shell_mass, -psi);
    }

    unit_velocity velocity(unit_mass m) const
    {
        auto u = gamma_beta(m);
        return u / std::sqrt(1.0 + u * u) * light_speed;
    }

    inverse<unit_mass> dudm(unit_mass m) const
    {
        return -psi / m * gamma_beta(m);
    }


    //=========================================================================
    unit_length radius(unit_mass m, unit_time t) const
    {
        return velocity(m) * t;
    }

    unit_mass_density density(unit_mass m, unit_time t) const
    {
        double gamma_squared = 1.0 + pow<2>(gamma_beta(m));
        double beta = velocity(m) / light_speed;
        return gamma_squared * beta / (solid_angle * pow<3>(radius(m, t))) / abs(dudm(m));
    }


    //=========================================================================
    unit_scalar cloud_gamma_beta() const
    {
        auto beta = velocity(envelop_mass) / light_speed;
        return beta / std::sqrt(1.0 - beta * beta);
    }

    unit_velocity cloud_velocity() const
    {
        return velocity(envelop_mass);
    }


    //=========================================================================
    unit_length cloud_outer_boundary(unit_time t) const
    {
        return cloud_velocity() * t;
    }

    unit_length envelop_outer_boundary(unit_time t) const
    {
        return radius(fastest_shell_mass, t);
    }

    unit_mass cloud_mass(unit_time t) const
    {
        double n1 = 2.0;
        auto r0 = inner_radius;
        auto rc = cloud_outer_boundary(t);
        return n1 == 3.0
        ? solid_angle * (density_at(rc, t) * pow<3>(rc) * std::log(rc / r0))
        : solid_angle * (density_at(rc, t) * pow<3>(rc) - density_at(r0, t) * pow<3>(r0)) / (3 - n1);
    }

    unit_mass total_mass(unit_time t) const
    {
        return cloud_mass(t) + envelop_mass;
    }


    //=========================================================================
    unit_mass mass_coordinate(unit_length r, unit_time t) const
    {
        auto f = [this, r, t] (unit_mass m)
        {
            return std::log10(r.value) - std::log10(radius(m, t).value);
        };
        return detail::solve_secant(f, fastest_shell_mass, fastest_shell_mass * 2, 1e-10);
    }

    unit_mass_density power_law_density_at(unit_length r, unit_time t) const
    {
        auto r_outer = cloud_outer_boundary(t);
        auto d_outer = density_at(r_outer, t);
        return d_outer * std::pow(r / r_outer, -2.0);
    }

    unit_mass_density density_at(unit_length r, unit_time t) const
    {
        auto r1 = envelop_outer_boundary(t);

        if (r < cloud_outer_boundary(t))
            return power_law_density_at(r, t);
        if (r > r1)
            return density_at(r1, t) * std::pow(r / r1, -2.0);
        return density(mass_coordinate(r, t), t);
    }

    unit_scalar gamma_beta_at(unit_length r, unit_time t) const
    {
        auto r1 = envelop_outer_boundary(t);

        if (r < cloud_outer_boundary(t))
            return cloud_gamma_beta();
        if (r > r1)
            return gamma_beta(mass_coordinate(r1, t));
        return gamma_beta(mass_coordinate(r, t));
    }

    unit_velocity velocity_at(unit_length r, unit_time t) const
    {
        auto u = gamma_beta_at(r, t);
        return u / std::sqrt(1.0 + u * u) * light_speed;
    }


    //=============================================================================
    auto primitive_srhd(unit_length r, unit_time t, unit_scalar entropy=-6.0) const
    {
        auto d = density_at(r, t);
        auto u = gamma_beta_at(r, t);
        auto p = unit_energy_density(1e-8);
        // auto p = unit_energy_density(std::pow(d.value, 4. / 3) * std::exp(entropy));
        return numeric::tuple(d, u, unit_scalar(0.0), unit_scalar(0.0), p);
    }


    //=========================================================================
    unit_length   inner_radius       = 1.0;
    unit_mass     envelop_mass       = 1.0;
    unit_mass     fastest_shell_mass = 1e-5;
    unit_scalar   fastest_gamma_beta = 4.0;
    unit_scalar   solid_angle        = 4 * M_PI;
    unit_velocity light_speed        = 1.0;
    double        psi                = 0.25;
};

} // namespace mara




//=============================================================================
#define DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_atmospheres()
{
}
