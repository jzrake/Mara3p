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




#include "core_ndarray.hpp"
#include "core_dimensional.hpp"
#include "mesh_sliding.hpp"
#include "physics_srhd.hpp"




//=============================================================================
namespace sedov {




//=============================================================================
struct radial_track_t
{
    nd::shared_array<dimensional::unit_length, 1> face_radii;
    dimensional::unit_scalar                      theta0;
    dimensional::unit_scalar                      theta1;
};

using primitive_per_length_t = decltype(srhd::primitive_t() / dimensional::unit_length());
using primitive_per_scalar_t = decltype(srhd::primitive_t() / dimensional::unit_scalar());
using radial_godunov_data_t  = std::tuple<srhd::flux_vector_t, dimensional::unit_velocity>;
using polar_godunov_data_t   = std::tuple<srhd::flux_vector_t, dimensional::unit_area, nd::uint, nd::uint>;
using track_data_t           = std::tuple<radial_track_t, nd::shared_array<srhd::primitive_t, 1>, nd::shared_array<primitive_per_length_t, 1>>;
using primitive_function_t   = std::function<srhd::primitive_t(dimensional::unit_length, dimensional::unit_scalar)>;




//=============================================================================
dimensional::unit_area face_area(
    dimensional::unit_length r0, dimensional::unit_length r1,
    dimensional::unit_scalar t0, dimensional::unit_scalar t1);

dimensional::unit_volume cell_volume(
    dimensional::unit_length r0, dimensional::unit_length r1,
    dimensional::unit_scalar t0, dimensional::unit_scalar t1);




//=============================================================================
dimensional::unit_length                      minimum_spacing  (radial_track_t track);
dimensional::unit_scalar                      cell_center_theta(radial_track_t track);
nd::shared_array<dimensional::unit_length, 1> cell_center_radii(radial_track_t track);
nd::shared_array<dimensional::unit_area,   1> radial_face_areas(radial_track_t track);
nd::shared_array<dimensional::unit_volume, 1> cell_volumes     (radial_track_t track);




/**
 * @brief      Convenience function to build a radial track with logarithmically
 *             spaced, reasonably isotropic cells.
 *
 * @param      r0      The inner radius
 * @param      r1      The outer radius
 * @param[in]  theta0  The lower polar angle
 * @param[in]  theta1  The upper polar angle
 *
 * @return     A radial track
 */
radial_track_t generate_radial_track(
    dimensional::unit_length r0,
    dimensional::unit_length r1,
    dimensional::unit_scalar theta0, 
    dimensional::unit_scalar theta1);




/**
 * @brief      Return an array of conserved quantities in the given track
 *
 * @param      track      The radial track itself (geometric data)
 * @param[in]  primitive  A mapping (r, theta) -> primitive
 *
 * @return     A shared array of conserved data
 */
nd::shared_array<srhd::conserved_t, 1> generate_conserved(radial_track_t track, primitive_function_t primitive);




/**
 * @brief      Recover the primtive data on a track from conserved quantities.
 *
 * @param      track  The radial track itself (geometric data)
 * @param      uc     The array of conserved quantities
 *
 * @return     A shared array of primitive data
 */
nd::shared_array<srhd::primitive_t, 1> recover_primitive(
    radial_track_t track,
    nd::shared_array<srhd::conserved_t, 1> uc);




/**
 * @brief      Compute PLM-estimated gradient of primitive quantities in the
 *             radial direction.
 *
 * @param      track  The radial track itself (geometric data)
 * @param      pc     The array of primitive data in the cells on that track
 *
 * @return     A shared array of primitive gradients
 *
 * @note       The gradients in the first and last cells are set to zero, and
 *             the array of gradients has the same shape as the array of
 *             primitive data.
 */
nd::shared_array<primitive_per_length_t, 1> radial_gradient(
    radial_track_t track,
    nd::shared_array<srhd::primitive_t, 1> pc);




/**
 * @brief      Generate Godunov data for the radially oriented faces, consisting
 *             of the radial HLLC flux across the contact discontinuity, and the
 *             radial velocity of the contact.
 *
 * @param      track  The radial track itself (geometric data)
 * @param      pc     The array of primitive data in the cells on that track
 * @param      dc     The radial derivative of the primitive quantities
 *
 * @return     A shared array of radial Godunov data
 */
nd::shared_array<radial_godunov_data_t, 1> radial_godunov_data(
        radial_track_t track,
        nd::shared_array<srhd::primitive_t, 1> pc,
        nd::shared_array<primitive_per_length_t, 1> dc,
        radial_godunov_data_t inner_boundary_data,
        radial_godunov_data_t outer_boundary_data);




/**
 * @brief      Generate Godunov data for the theta-oriented faces between two
 *             radial tracks, consisting of the area-weighted flux and the
 *             indexes (in their respective radial tracks) of the cells to
 *             either side of the face.
 *
 * @param[in]  t0    The track two to the left
 * @param[in]  t1    The track to the left of the interface
 * @param[in]  t2    The track to the right of the interface
 * @param[in]  t3    The track two to the right
 *
 * @return     A shared array of polar Godunov data
 *
 * @note       If the data in either of t0 or t3 is empty, then only t1 and t2
 *             are used, and piecewise constant extrapolation is used in the
 *             polar direction.
 */
nd::shared_array<polar_godunov_data_t, 1> polar_godunov_data(track_data_t t0, track_data_t t1, track_data_t t2, track_data_t t3);




/**
 * @brief      Return the difference of conserved quantities in the cells of a
 *             track, given pre-computed Godunov data and the length of the time
 *             step.
 *
 * @param      track  The radial track itself (geometric data)
 * @param      ff     The radial Godunov data on the track
 * @param      gfl    The polar Godunov data on the track's left boundary
 * @param      gfr    The polar Godunov data on the track's right boundary
 * @param      dt     The time step
 *
 * @return     A shared array of conserved quantity differences
 */
nd::shared_array<srhd::conserved_t, 1> delta_conserved(radial_track_t track,
    nd::shared_array<srhd::primitive_t, 1> pc,
    nd::shared_array<radial_godunov_data_t, 1> ff,
    nd::shared_array<polar_godunov_data_t, 1> gfl,
    nd::shared_array<polar_godunov_data_t, 1> gfr,
    dimensional::unit_time dt);




/**
 * @brief      Return the difference in radial face positions given pre-computed
 *             Godunov data and the length of the time step.
 *
 * @param      ff    The radial Godunov data on the track
 * @param      dt    The time step
 *
 * @return     A shared array of radial positions
 */
nd::shared_array<dimensional::unit_length, 1> delta_face_positions(
    nd::shared_array<radial_godunov_data_t, 1> ff,
    dimensional::unit_time dt);




/**
 * @brief      Sample the primitive quantities at a given radius on a track,
 *             using full track data including the radial gradient of primitive
 *             quantities. The given index is used as a hint for where to begin
 *             searching for the cell that contains the target radius.
 *
 * @param[in]  track_data  The full track data
 * @param[in]  r           The radius at which to sample the track
 * @param[in]  index       The hint index
 * @param[in]  fallback    The fallback primitive state if the sample cannot be
 *                         made
 *
 * @return     An extrapolated primitive data
 *
 * @note       If the target radius lies either above or below the track bounds,
 *             then the fallback value is used instead.
 */
srhd::primitive_t sample(track_data_t track_data, dimensional::unit_length r, nd::uint index, srhd::primitive_t fallback);




/**
 * @brief      Extend a radial track and its primitive quantities by one zone in
 *             the radial direction.
 *
 * @param[in]  track                      The track
 * @param[in]  uc                         { parameter_description }
 * @param[in]  maximum_cell_aspect_ratio  The maximum cell aspect ratio
 *
 * @return     A pair of extended track and conserved data
 */
std::pair<radial_track_t, nd::shared_array<srhd::conserved_t, 1>> refine(
    radial_track_t track,
    nd::shared_array<srhd::conserved_t, 1> uc,
    dimensional::unit_scalar maximum_cell_aspect_ratio);

} // namespace sedov
