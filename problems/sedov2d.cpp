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




#include "core_util.hpp"
#include "scheme_sedov2d.hpp"




//=============================================================================
int main()
{
    auto num_tracks = 100;

    // initial condition:
    auto t0 = sedov::generate_radial_tracks(num_tracks, 1.0, 10.0);
    auto u0 = t0 | nd::map(sedov::generate_conserved);

    // time step:
    auto p0 = nd::zip(t0, u0) | nd::map(util::apply_to(sedov::recover_primitive));
    auto [te, pc] = nd::unzip(nd::zip(t0, p0) | nd::map(util::apply_to(sedov::extend)));

    auto dc = nd::zip(te, pc) | nd::map(util::apply_to(sedov::radial_gradient));
    auto ff = nd::zip(te, pc, dc) | nd::map(util::apply_to(sedov::radial_godunov_data));
    auto gf = nd::zip(te, pc, dc) | nd::adjacent_zip() | nd::map(util::apply_to(sedov::polar_godunov_data));

    auto dt = dimensional::unit_time(1.0);
    auto du = sedov::delta_conserved(t0(0), p0(0), ff(0), gf(0), gf(1), dt);
    auto dr = sedov::delta_face_positions(ff(0), dt);


    u0(0) + du;
    t0(0).face_radii + dr;


    // auto [vertices, indexes] = sedov::quad_mesh(t0);
    // vtk::write(std::cout, "Grid", vertices, indexes, std::pair("cell radius", rc));


    return 0;
}
