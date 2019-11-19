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
#include "app_state_templates.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "scheme_plm_gradient.hpp"




namespace mara {




//=============================================================================
template<
    typename ConservedType,
    typename PrimitiveType,
    typename RiemannSolverType,
    typename RecoverPrimitiveType,
    typename SourceTermsType,
    typename MeshGeometryType>

state_with_vertices_t<ConservedType> advance(

    state_with_vertices_t<ConservedType>       state,
    dimensional::unit_time                     dt,
    PrimitiveType                              inner_boundary_primitive,
    PrimitiveType                              outer_boundary_primitive,
    RiemannSolverType                          riemann_solver,
    RecoverPrimitiveType                       recover_primitive,
    SourceTermsType                            source_terms,
    MeshGeometryType                           mesh_geometry,
    double                                     plm_theta)
{
    try {
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

        auto [ibf, ibv] = riemann_solver(std::pair(inner_boundary_primitive, front(p0)));
        auto [obf, obv] = riemann_solver(std::pair(back(p0), outer_boundary_primitive));
        auto pl = select(p0 + 0.5 * dx * gx, 0, 0, -1);
        auto pr = select(p0 - 0.5 * dx * gx, 0, 1);
        auto fv = nd::zip(pl, pr)
        | nd::map(riemann_solver)
        | nd::extend_uniform_lower(std::pair(ibf, ibv))
        | nd::extend_zero_gradient_upper() //nd::extend_uniform_upper(std::pair(obf, obv))
        | nd::to_shared();

        auto ff = fv | nd::map([] (auto a) { return std::get<0>(a); });
        auto vf = fv | nd::map([] (auto a) { return std::get<1>(a); });
        auto df = ff | nd::multiply(da) | nd::adjacent_diff() | nd::to_shared();
        auto x1 = state.vertices  + (vf * dt)                 | nd::to_shared();
        auto u1 = state.conserved + (s0 - df) * dt            | nd::to_shared();

        return {
            state.iteration + 1,
            state.time + dt,
            x1,
            u1,
        };
    }
    catch (const std::exception& e)
    {
        auto x0 = mesh_geometry.cell_centers(state.vertices);
        auto dv = mesh_geometry.cell_volumes(state.vertices);

        for (auto [x, u] : nd::zip(x0, state.conserved / dv))
        {
            try {
                recover_primitive(u);
            }
            catch (const std::exception& f)
            {
                throw std::invalid_argument(std::string(f.what()) + " at r = " + std::to_string(x.value));
            }
        }
        throw;
    }
}

} // namespace mara
