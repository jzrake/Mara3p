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
struct cold_relativistic_wind_t
{


    //=========================================================================
    using mass_loss_rate_function_t = std::function<unit_mass_rate(unit_time)>;
    using gamma_beta_function_t     = std::function<unit_scalar(unit_time)>;


    //=========================================================================
    auto with_solid_angle(unit_scalar new_solid_angle) const
    -> cold_relativistic_wind_t
    {
        return {new_solid_angle, light_speed, mass_loss_rate, gamma_beta};
    }

    auto with_mass_loss_rate(mass_loss_rate_function_t new_mass_loss_rate) const
    -> cold_relativistic_wind_t
    {
        return {solid_angle, light_speed, new_mass_loss_rate, gamma_beta};
    }

    auto with_gamma_beta(gamma_beta_function_t new_gamma_beta) const
    -> cold_relativistic_wind_t
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

    auto primitive(unit_length r, unit_time t, unit_scalar temperature=1e-3) const
    {
        auto d = density_at(r, t);
        auto u = gamma_beta(t);
        auto p = d * pow<2>(light_speed) * temperature;
        return numeric::tuple(d, u, unit_scalar(0.0), unit_scalar(0.0), p);
    }


    //=========================================================================
    unit_scalar   solid_angle = 4 * M_PI;
    unit_velocity light_speed = 1.0;
    mass_loss_rate_function_t mass_loss_rate;
    gamma_beta_function_t     gamma_beta;
};




/**
 * @brief      A helper function that models a fast-rising pulse, with amplitude
 *             of 1.0 and a given duration and onset time.
 *
 * @param[in]  onset_time  The onset time
 * @param[in]  decay_time  The decay time
 *
 * @return     A function modeling a FRED-like pulse profile
 */
inline auto fast_rise_exponential_decay(dimensional::unit_time onset_time, dimensional::unit_time decay_time)
{
    return [=] (dimensional::unit_time t)
    {
        auto x = (t - onset_time) / decay_time;
        return x < 0.0 ? 0.0 : 0.25 * std::pow(x, 2.0) * std::exp(2.0 - x);
        // return x < 0.0 ? 0.0 : std::sqrt(2.0 * x) * std::exp(0.5 - x);
    };
}

} // namespace mara




//=============================================================================
#define DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_atmospheres()
{
}
