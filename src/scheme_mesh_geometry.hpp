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
#include "core_dimensional.hpp"
#include "core_ndarray_ops.hpp"
#include "core_util.hpp"




//=============================================================================
namespace mara {




//=============================================================================
struct spherical_mesh_geometry_t
{
    using vertices_array_t = nd::shared_array<dimensional::unit_length, 1>;

    static auto shell_volume(dimensional::unit_length r0, dimensional::unit_length r1)
    {
        return (r1 * r1 * r1 - r0 * r0 * r0) / 3.0;
    }

    static auto cell_centers(vertices_array_t vertices)
    {
        return vertices | nd::adjacent_mean();
    }

    static auto cell_spacing(vertices_array_t vertices)
    {
        return vertices | nd::adjacent_diff();
    }

    static auto cell_volumes(vertices_array_t vertices)
    {
        return vertices | nd::adjacent_zip() | nd::map(util::apply_to(shell_volume));
    }

    static auto face_areas(vertices_array_t vertices)
    {
        return vertices | nd::map([] (auto r) { return r * r; });
    }
};

} // namespace mara
