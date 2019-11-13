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
#include "core_ndarray_ops.hpp"
#include "scheme_plm_gradient.hpp"




namespace mara {




//=============================================================================
template<
    typename ConservedType,
    typename FluxVectorType,
    typename RiemannSolverType,
    typename RecoverPrimitiveType,
    typename SourceTermsType,
    typename MeshGeometryType>

state_with_vertices_t<ConservedType> advance(

    state_with_vertices_t<ConservedType>       state,
    dimensional::unit_time                     dt,
    FluxVectorType                             inner_boundary_flux,
    RiemannSolverType                          riemann_solver,
    RecoverPrimitiveType                       recover_primitive,
    SourceTermsType                            source_terms,
    MeshGeometryType                           mesh_geometry,
    double                                     plm_theta)
{
    auto bf = inner_boundary_flux;
    auto bv = dimensional::unit_velocity(0.0);
    auto da = mesh_geometry.face_areas(state.vertices);
    auto dv = mesh_geometry.cell_volumes(state.vertices);
    auto dx = mesh_geometry.cell_spacing(state.vertices);
    auto x0 = mesh_geometry.cell_centers(state.vertices);
    auto p0 = state.conserved / dv | nd::map(recover_primitive) | nd::to_shared();
    auto s0 = nd::zip(p0, x0) | nd::map(source_terms) | nd::multiply(dv);
    auto gx = nd::zip(x0 | nd::adjacent_zip3(), p0 | nd::adjacent_zip3())
    | nd::map(mara::plm_gradient(plm_theta))
    | nd::extend_zeros()
    | nd::to_shared();

    auto pl = select(p0 + 0.5 * dx * gx, 0, 0, -1);
    auto pr = select(p0 - 0.5 * dx * gx, 0, 1);
    auto fv = nd::zip(pl, pr)
    | nd::map(riemann_solver)
    | nd::extend_uniform_lower(std::pair(bf, bv))
    | nd::extend_zero_gradient_upper()
    | nd::to_shared();

    auto ff = fv | nd::map([] (auto a) { return std::get<0>(a); });
    auto vf = fv | nd::map([] (auto a) { return std::get<1>(a); });
    auto df = ff | nd::multiply(da) | nd::adjacent_diff() | nd::to_shared();
    auto q1 = state.conserved + (s0 - df) * dt | nd::to_shared();
    auto x1 = state.vertices  + (vf * dt) | nd::to_shared();

    return {
        state.iteration + 1,
        state.time + dt,
        x1,
        q1,
    };
}

} // namespace mara
