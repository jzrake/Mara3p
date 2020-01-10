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
#include "core_bqo_tree.hpp"
#include "core_geometric.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_rational.hpp"
#include "core_sequence.hpp"
#include "mesh_cartesian_3d.hpp"
#include "physics_mhd.hpp"




//=============================================================================
namespace mhd_scheme_v2 {




//=============================================================================
using position_t                   = geometric::euclidean_vector_t<dimensional::unit_length>;
using multilevel_index_t           = bqo_tree::tree_index_t<3>;
using cell_primitive_variables_t   = nd::shared_array<mhd::primitive_t, 3>;
using cell_conserved_density_t     = nd::shared_array<mhd::conserved_density_t, 3>;
using face_godunov_data_t          = std::array<nd::shared_array<mhd::godunov_data_t, 3>, 3>;
using face_magnetic_flux_density_t = std::array<nd::shared_array<mhd::unit_magnetic_field, 3>, 3>;
using edge_electromotive_density_t = std::array<nd::shared_array<mhd::unit_electric_field, 3>, 3>;
using vector_potential_function_t  = std::function<mhd::vector_potential_t(position_t)>;
using primitive_function_t         = std::function<mhd::primitive_t(position_t, mhd::magnetic_field_vector_t)>;




/**
 * @brief      Recover an array of primitive variables from the cell-like
 *             conserved densities and face-like magnetic flux densities.
 *
 * @param[in]  uc    An array of conserved densities
 * @param[in]  bf    A 3-array of face-like magnetic flux densities
 *
 * @return     An array of primitive variables
 */
cell_primitive_variables_t primitive_array(cell_conserved_density_t uc, face_magnetic_flux_density_t bf);




/**
 * @brief      Construct an array of edge-like, vector-potential data, for the
 *             given multi-level block index. Vertices are assumed to be
 *             arranged to cover the unit lattice on level 0.
 *
 * @param[in]  index             The block index
 * @param[in]  block_size        The block size (isotropic)
 * @param[in]  vector_potential  The vector potential function
 *
 * @return     A 3-array of edge-like "vector potential" data
 */
edge_electromotive_density_t construct_vector_potential(multilevel_index_t index, unsigned block_size, vector_potential_function_t vector_potential);




/**
 * @brief      Construct an array of edge-like, vector-potential data, for the
 *             given multi-level block index. Vertices are assumed to be
 *             arranged to cover the unit lattice on level 0.
 *
 * @param[in]  index               The block index
 * @param[in]  block_size          The block size (isotropic)
 * @param[in]  bf                  A 3-array of face-like magnetic flux
 *                                 densities
 * @param[in]  primitive_function  The primitive variable function
 *
 * @return     An array of cell-like conserved densities
 */
cell_conserved_density_t construct_conserved(multilevel_index_t index, unsigned block_size, face_magnetic_flux_density_t bf, primitive_function_t primitive_function);




/**
 * @brief      Construct an extended array of cell-like primitive quantities.
 *
 * @param[in]  pc_vector   A vector of cell-like primitive variables
 * @param[in]  block_size  The block size of each of the arrays
 * @param[in]  count       The number of zones to add on  all axes
 *
 * @return     An extended primitive variable array
 */
cell_primitive_variables_t extend_cell_primitive_variables(std::vector<cell_primitive_variables_t> pc_vector, nd::uint block_size, nd::uint count);




/**
 * @brief      Construct an extended array of face-like magnetic flux densities.
 *
 * @param[in]  bf_vector   A vector of face-like magnetic field arrays
 * @param[in]  block_size  The block size of the associated finite volume mesh
 * @param[in]  count       The number of zones to add on the transverse axes
 *
 * @return     A 3-array of extended face-centered magnetic flux densities
 */
face_magnetic_flux_density_t extend_face_magnetic_flux_density(std::vector<face_magnetic_flux_density_t> bf_vector, nd::uint block_size, nd::uint count);




/**
 * @brief      Generate Godunov fluxes and EMF's using piecewise constant
 *             (first-order) extrapolation.
 *
 * @param[in]  pc    Cell-centered primitive quantities, extended by 1 zones on
 *                   all three axes
 * @param[in]  bf    Face B-fields, extended by 1 zone on the transverse axes
 *
 * @return     A 3-array of face-centered Godunov fluxes and EMF's. Each of the
 *             returned arrays are extended by 1 zone on the transverse axes
 */
face_godunov_data_t godunov_fluxes(cell_primitive_variables_t pc, face_magnetic_flux_density_t bf);




/**
 * @brief      Generate EMF values on cell edges by averaging the electric field
 *             estimates from godunov data on the cell faces.
 *
 * @param[in]  gf    The godunoc data on cell faces
 *
 * @return     A 3-array of electric field values at edge midpoints, in the
 *             direction of the respective edge
 */
edge_electromotive_density_t electromotive_forces(face_godunov_data_t gf);




/**
 * @brief      Generate a new conserved density array from Godunov fluxes, for a
 *             given time step dt and uniform, isotropic mesh spacing dl.
 *
 * @param[in]  uc    The array of conserved densities
 * @param[in]  gf    The Godunov fluxes
 * @param[in]  dt    The time step
 * @param[in]  dl    The mesh spacing
 *
 * @return     A new array of conserved densities
 */
cell_conserved_density_t updated_conserved_density(cell_conserved_density_t uc, face_godunov_data_t gf, dimensional::unit_time dt, dimensional::unit_length dl);




/**
 * @brief      Evaluate magnetic fields on cell faces by computing the staggered
 *             curl of electric field data on cell edges, as if advancing
 *             Faraday's law by a unit length of time.
 *
 * @param[in]  ee    The electric field on cell edges
 * @param[in]  dl    The uniform, isotropic mesh spacing
 *
 * @return     A 3-array of extended face-centered magnetic flux densities
 */
face_magnetic_flux_density_t curl(edge_electromotive_density_t ee, dimensional::unit_length dl);




/**
 * @brief      Generate a new magnetic flux density array from edge EMF's, for a
 *             given time step dt and uniform, isotropic mesh spacing dl.
 *
 * @param[in]  bf    A 3-array of face-like magnetic flux densities
 * @param[in]  ee    A 3-array of edge-like electromotive force densities
 * @param[in]  dt    The time step
 * @param[in]  dl    The mesh spacing
 *
 * @return     A new array of magnetic flux densities
 */
face_magnetic_flux_density_t updated_magnetic_flux_density(face_magnetic_flux_density_t bf, edge_electromotive_density_t ee, dimensional::unit_time dt, dimensional::unit_length dl);

} // namespace mara
