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




//=============================================================================
dimensional::unit_area face_area(
    dimensional::unit_length r0, dimensional::unit_length r1,
    dimensional::unit_scalar t0, dimensional::unit_scalar t1);

dimensional::unit_volume cell_volume(
    dimensional::unit_length r0, dimensional::unit_length r1,
    dimensional::unit_scalar t0, dimensional::unit_scalar t1);




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
 * @param      track  The radial track itself (geometric data)
 *
 * @return     A shared array of conserved data
 */
nd::shared_array<srhd::conserved_t, 1> generate_conserved(radial_track_t track);




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
        nd::shared_array<primitive_per_length_t, 1> dc);




/**
 * @brief      Generate Godunov data for the theta-oriented faces between two
 *             radial tracks, consisting of the area-weighted flux and the
 *             indexes (in their respective radial tracks) of the cells to
 *             either side of the face.
 *
 * @param      L     The track on the left
 * @param      R     The track on the right
 *
 * @return     A shared array of polar Godunov data
 */
nd::shared_array<polar_godunov_data_t, 1> polar_godunov_data(track_data_t L, track_data_t R);




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
 * @brief      Extend a radial track and its primitive quantities by one zone in
 *             the radial direction.
 *
 * @param      track  The track to extend
 * @param      pc     The primitive data
 *
 * @return     A pair of extended track and primitive data
 */
std::pair<radial_track_t, nd::shared_array<srhd::primitive_t, 1>> extend(
    radial_track_t track,
    nd::shared_array<srhd::primitive_t, 1> pc);




//=============================================================================
inline auto radial_face_areas(radial_track_t track)
{
    return track.face_radii
    | nd::map([t0=track.theta0, t1=track.theta1] (auto r)
    {
        return face_area(r, r, t0, t1);
    });
}

inline auto cell_volumes(radial_track_t track)
{
    return track.face_radii
    | nd::adjacent_zip()
    | nd::map(util::apply_to([t0=track.theta0, t1=track.theta1] (auto r0, auto r1)
    {
        return cell_volume(r0, r1, t0, t1);
    }));
}

inline auto cell_center_radii(radial_track_t track)
{
    return track.face_radii | nd::adjacent_mean();
}

inline auto cell_center_theta(radial_track_t track)
{
    return 0.5 * (track.theta0 + track.theta1);
}

inline auto polar_faces(radial_track_t L, radial_track_t R)
{
    return mesh::transverse_faces(L.face_radii, R.face_radii);
}

inline auto minimum_spacing(radial_track_t track)
{
    return nd::min(track.face_radii | nd::adjacent_diff());
}

} // namespace sedov
