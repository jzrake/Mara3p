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
#include <algorithm>
#include "core_dimensional.hpp"
#include "core_dimensional_math.hpp"
#include "core_geometric.hpp"
#include "core_numeric_tuple.hpp"




//=============================================================================
namespace two_body {


//=============================================================================
using namespace dimensional;




//=============================================================================
struct point_mass_t
{
    unit_mass     mass       = 1.0;
    unit_length   position_x = 0.0;
    unit_length   position_y = 0.0;
    unit_velocity velocity_x = 0.0;
    unit_velocity velocity_y = 0.0;
};




//=============================================================================
using orbital_state_t       = std::pair<point_mass_t, point_mass_t>;
using orbital_elements_t    = numeric::tuple_t<unit_length, unit_mass, unit_scalar, unit_scalar>;
using orbital_orientation_t = numeric::tuple_t<unit_length, unit_length, unit_velocity, unit_velocity, unit_scalar, unit_time>;
using orbital_parameters_t  = numeric::tuple_t<orbital_elements_t, orbital_orientation_t>;
const static auto G         = unit_energy(1.0) * unit_length(1.0) / unit_mass(1.0) / unit_mass(1.0);




//=============================================================================
inline auto semimajor_axis(orbital_elements_t el) { return numeric::get<0>(el); }
inline auto total_mass    (orbital_elements_t el) { return numeric::get<1>(el); }
inline auto mass_ratio    (orbital_elements_t el) { return numeric::get<2>(el); }
inline auto eccentricity  (orbital_elements_t el) { return numeric::get<3>(el); }

inline auto orbital_elements(unit_length a=1.0, unit_mass M=1.0, unit_scalar q=1.0, unit_scalar e=0.0)
{
    return orbital_elements_t{{a, M, q, e}};
}

inline auto omega(orbital_elements_t el)
{
    auto M = total_mass(el);
    auto a = semimajor_axis(el);
    return sqrt(G * M / a / a / a);
}

inline auto period(orbital_elements_t el)
{
    return 2 * M_PI / omega(el);
}

inline auto eccentric_anomaly(orbital_elements_t el, unit_time t)
{
    auto solve_newton_rapheson = [] (auto f, auto g, double x)
    {
        while (std::abs(f(x)) > 1e-12)
        {
            x -= f(x) / g(x);
        }
        return x;
    };

    auto e = eccentricity(el);
    auto M = double(omega(el) * t);
    auto f = [e, M] (double E) { return E - e * std::sin(E) - M; };
    auto g = [e]    (double E) { return 1 - e * std::cos(E); };

    return solve_newton_rapheson(f, g, M);
}

inline auto orbital_state_from_eccentric_anomaly(orbital_elements_t el, unit_scalar eccentric_anomaly)
{
    auto a  = semimajor_axis(el);
    auto q  = mass_ratio(el);
    auto e  = eccentricity(el);
    auto W  = omega(el);
    auto mu = 1.0 / (1.0 + q);
    auto M1 = total_mass(el) * (1.0 - mu);
    auto M2 = total_mass(el) * mu;
    auto _E = eccentric_anomaly;

    auto x1  = -a * mu * (e - std::cos(_E));
    auto y1  = +a * mu * (0 + std::sin(_E)) * std::sqrt(1 - e * e);
    auto x2  = -x1 / q;
    auto y2  = -y1 / q;
    auto vx1 = -a * mu * W / (1 - e * std::cos(_E)) * std::sin(_E);
    auto vy1 = +a * mu * W / (1 - e * std::cos(_E)) * std::cos(_E) * std::sqrt(1 - e * e);
    auto vx2 = -vx1 / q;
    auto vy2 = -vy1 / q;

    auto c1 = point_mass_t{M1, x1, y1, vx1, vy1};
    auto c2 = point_mass_t{M2, x2, y2, vx2, vy2};

    return orbital_state_t{c1, c2};
}

inline auto orbital_state(orbital_elements_t el, unit_time t)
{
    return orbital_state_from_eccentric_anomaly(el, eccentric_anomaly(el, t));
}




//=============================================================================
inline auto cm_position_x    (orbital_orientation_t o) { return numeric::get<0>(o); }
inline auto cm_position_y    (orbital_orientation_t o) { return numeric::get<1>(o); }
inline auto cm_velocity_x    (orbital_orientation_t o) { return numeric::get<2>(o); }
inline auto cm_velocity_y    (orbital_orientation_t o) { return numeric::get<3>(o); }
inline auto periapse_argument(orbital_orientation_t o) { return numeric::get<4>(o); }
inline auto periapse_time    (orbital_orientation_t o) { return numeric::get<5>(o); }




//=============================================================================
inline auto kinetic_energy(point_mass_t p)
{
    auto vx = p.velocity_x;
    auto vy = p.velocity_y;
    return 0.5 * p.mass * (vx * vx + vy * vy);
}

inline auto potential(point_mass_t p, unit_length x, unit_length y, unit_length softening_length)
{
    auto dx = p.position_x - x;
    auto dy = p.position_y - y;
    auto r2 = dx * dx + dy * dy;
    auto s2 = softening_length * softening_length;
    return -G * p.mass / sqrt(r2 + s2);
}

inline auto perturb(point_mass_t p, unit_mass dM)
{
    p.mass += dM;
    return p;
}

inline auto perturb(point_mass_t p, unit_momentum dpx, unit_momentum dpy)
{
    p.velocity_x += dpx / p.mass;
    p.velocity_y += dpy / p.mass;
    return p;
}

inline auto perturb(point_mass_t p, unit_mass dM, unit_momentum dpx, unit_momentum dpy)
{
    p.velocity_x += dpx / p.mass;
    p.velocity_y += dpy / p.mass;
    p.mass += dM;
    return p;
}




//=============================================================================
inline auto total_mass(orbital_state_t s)
{
    return s.first.mass + s.second.mass;
}

inline auto mass_ratio(orbital_state_t s)
{
    return s.second.mass / s.first.mass;
}

inline auto separation(orbital_state_t s)
{
    auto x1 = s.first.position_x;
    auto y1 = s.first.position_y;
    auto x2 = s.second.position_x;
    auto y2 = s.second.position_y;
    return sqrt(pow<2>(x2 - x1) + pow<2>(y2 - y1));
}

inline auto kinetic_energy(orbital_state_t s)
{
    return kinetic_energy(s.first) + kinetic_energy(s.second);
}

inline auto potential(orbital_state_t s, unit_length x, unit_length y, unit_length softening_length)
{
    return potential(s.first, x, y, softening_length) + potential(s.second, x, y, softening_length);
}

inline auto total_energy(orbital_state_t s)
{
    return kinetic_energy(s) - G * s.first.mass * s.second.mass / separation(s);
}

inline auto recover_orbital_elements(orbital_state_t state, unit_time t)
{
    auto [c1, c2] = state;


    // component masses, total mass, and mass ratio
    auto M1 = c1.mass;
    auto M2 = c2.mass;
    auto M = M1 + M2;
    auto q = M2 / M1;


    // position and velocity of the CM frame
    auto x_cm  = (c1.position_x * c1.mass + c2.position_x * c2.mass) / M;
    auto y_cm  = (c1.position_y * c1.mass + c2.position_y * c2.mass) / M;
    auto vx_cm = (c1.velocity_x * c1.mass + c2.velocity_x * c2.mass) / M;
    auto vy_cm = (c1.velocity_y * c1.mass + c2.velocity_y * c2.mass) / M;


    // positions and velocities of the components in the CM frame
    auto x1 = c1.position_x - x_cm;
    auto y1 = c1.position_y - y_cm;
    auto x2 = c2.position_x - x_cm;
    auto y2 = c2.position_y - y_cm;
    auto r1 = sqrt(x1 * x1 + y1 * y1);
    auto r2 = sqrt(x2 * x2 + y2 * y2);
    auto vx1 = c1.velocity_x - vx_cm;
    auto vy1 = c1.velocity_y - vy_cm;
    auto vx2 = c2.velocity_x - vx_cm;
    auto vy2 = c2.velocity_y - vy_cm;
    auto vf1 = -vx1 * y1 / r1 + vy1 * x1 / r1;
    auto vf2 = -vx2 * y2 / r2 + vy2 * x2 / r2;
    auto v1 = sqrt(vx1 * vx1 + vy1 * vy1);


    // energy and angular momentum
    auto E1 = 0.5 * M1 * (vx1 * vx1 + vy1 * vy1);
    auto E2 = 0.5 * M2 * (vx2 * vx2 + vy2 * vy2);
    auto L1 = M1 * r1 * vf1;
    auto L2 = M2 * r2 * vf2;
    auto R = r1 + r2;
    auto L = L1 + L2;
    auto E = E1 + E2 - G * M1 * M2 / R;

    if (E >= unit_energy(0.0))
    {
        throw std::invalid_argument("two_body::recover_orbital_elements (not a bound orbit)");
    }


    // Semi-major, semi-minor axes; eccentricity, apsides
    auto a = -0.5 * G * M1 * M2 / E;
    auto b = sqrt(-0.5 * L * L / E * (M1 + M2) / (M1 * M2));
    auto e = sqrt(std::clamp(1.0 - b * b / a / a, 0.0, 1.0));
    auto omega = sqrt(G * M / a / a / a);


    // semi-major and semi-minor axes of the primary
    auto a1 = a * q / (1.0 + q);
    auto b1 = b * q / (1.0 + q);


    // cos of nu and f: phase angle and true anomaly
    auto cn = e == 0.0 ? x1 / r1 : unit_scalar(1.0 - r1 / a1) / e;
    auto cf = a1 / r1 * (cn - e);


    // sin of nu and f
    auto sn = e == 0.0 ? y1 / r1 : (vx1 * x1 + vy1 * y1) / (e * v1 * r1) * std::sqrt(1.0 - e * e * cn * cn);
    auto sf = (b1 / r1) * sn;


    // cos and sin of eccentric anomaly
    auto cE = (e + cf)                    / (1.0 + e * cf);
    auto sE = std::sqrt(1.0 - e * e) * sf / (1.0 + e * cf);


    // mean anomaly and tau
    auto _E = std::atan2(sE, cE);
    auto _M = _E - e * sE;
    auto tau = t - _M / omega;


    // cartesian components of semi-major axis, and the argument of periapse
    auto ax = +(cn - e) * x1 + sn * std::sqrt(1.0 - e * e) * y1;
    auto ay = +(cn - e) * y1 - sn * std::sqrt(1.0 - e * e) * x1;
    auto pomega = atan2(ay, ax);


    // final result
    auto elements    = orbital_elements_t{{a, M, q, e}};
    auto orientation = orbital_orientation_t{{x_cm, y_cm, vx_cm, vy_cm, pomega, tau}};


    return orbital_parameters_t{{elements, orientation}};
}

inline auto transform(orbital_state_t state, orbital_orientation_t o)
{
    auto M1  = state.first.mass;
    auto x1  = state.first.position_x;
    auto y1  = state.first.position_y;
    auto vx1 = state.first.velocity_x;
    auto vy1 = state.first.velocity_y;
    auto M2  = state.second.mass;
    auto x2  = state.second.position_x;
    auto y2  = state.second.position_y;
    auto vx2 = state.second.velocity_x;
    auto vy2 = state.second.velocity_y;

    auto c = std::cos(-periapse_argument(o));
    auto s = std::sin(-periapse_argument(o));

    auto x1p  = +x1  * c + y1  * s + cm_position_x(o);
    auto y1p  = -x1  * s + y1  * c + cm_position_y(o);
    auto x2p  = +x2  * c + y2  * s + cm_position_x(o);
    auto y2p  = -x2  * s + y2  * c + cm_position_y(o);
    auto vx1p = +vx1 * c + vy1 * s + cm_velocity_x(o);
    auto vy1p = -vx1 * s + vy1 * c + cm_velocity_y(o);
    auto vx2p = +vx2 * c + vy2 * s + cm_velocity_x(o);
    auto vy2p = -vx2 * s + vy2 * c + cm_velocity_y(o);

    auto c1 = point_mass_t{M1, x1p, y1p, vx1p, vy1p};
    auto c2 = point_mass_t{M2, x2p, y2p, vx2p, vy2p};

    return orbital_state_t{c1, c2};
}

inline auto orbital_state(orbital_elements_t el, orbital_orientation_t o, unit_time t)
{
    return transform(orbital_state(el, t - periapse_time(o)), o);
}

} // namespace two_body




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
inline void test_two_body()
{
    using namespace two_body;

    auto t        = dimensional::unit_time(0.5);
    auto elements = orbital_elements(1.0, 1.0, 1.0, 0.5);
    auto state    = orbital_state(elements, t);

    require(abs(eccentricity  (numeric::get<0>(recover_orbital_elements(state, t))) - eccentricity  (elements)) < unit_scalar(1e-10));
    require(abs(semimajor_axis(numeric::get<0>(recover_orbital_elements(state, t))) - semimajor_axis(elements)) < unit_length(1e-10));
}

#endif // DO_UNIT_TESTS
