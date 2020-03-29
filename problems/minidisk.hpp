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
#include "app_config.hpp"
#include "core_bsp_tree.hpp"
#include "core_bqo_tree.hpp"
#include "core_ndarray.hpp"
#include "physics_iso2d.hpp"
#include "parallel_computable.hpp"




//=============================================================================
namespace minidisk {




//=============================================================================
using namespace dimensional;
using conserved_array_t = nd::shared_array<iso2d::conserved_density_t, 2>;
using primitive_array_t = nd::shared_array<iso2d::primitive_t, 2>;
using godunov_f_array_t = nd::shared_array<iso2d::flux_vector_t, 2>;
using conserved_tree_t  = bsp::shared_tree<conserved_array_t, 4>;
using unit_viscosity    = dimensional::quantity_t<0, 2,-1>;




//=============================================================================
struct side_effect_t
{
    unit_time next_due = 0.0;
    int count = 0;
};

struct schedule_t
{
    side_effect_t checkpoint;
    side_effect_t time_series;
};

struct solution_t
{
    rational::number_t                                      iteration;
    dimensional::unit_time                                  time;
    bsp::shared_tree<mpr::computable<conserved_array_t>, 4> conserved;
};

struct solver_data_t
{
    int                 block_size;
    int                 depth;
    unit_rate           omega_frame;
    unit_scalar         mach_number;
    unit_viscosity      kinematic_viscosity;
    unit_length         domain_radius;
    unit_length         softening_length;
    unit_scalar         eccentricity;
    unit_scalar         mass_ratio;
    unit_length         sink_radius;
    unit_rate           sink_rate;
    unit_rate           buffer_rate;
    unit_length         buffer_scale;
    unit_scalar         cfl;
};




//=============================================================================
inline auto config_template()
{
    return mara::config_template()
    .item("block_size",            64, "number of zones per side")
    .item("depth",                  1, "number of levels in the mesh")
    .item("omega_frame",          1.0, "reference frame rotation frequency (1.0 for co-rotating)")
    .item("mach_number",         10.0, "Mach number of the disk")
    .item("nu",                   0.0, "Kinematic viscosity coefficient")
    .item("domain_radius",        1.5, "half-size of the domain")
    .item("eccentricity",         0.0, "binary orbital eccentricity")
    .item("mass_ratio",           1.0, "binary mass ratio M2 / M1")
    .item("softening_length",    0.01, "gravitational softening length")
    .item("sink_radius",         0.01, "radius of mass sinks")
    .item("sink_rate",            1e3, "rate of mass and momentum removal at sinks")
    .item("buffer_rate",          1e2, "maximum rate of buffer driving")
    .item("buffer_scale",         0.2, "length scale of buffer onset")
    .item("cfl",                  0.3, "Courant number")
    .item("orbits",             100.0, "time to stop the simulation")
    .item("cpi",                  0.1, "number of orbits between checkpoints")
    .item("threads",                1, "number of threads to execute on")
    .item("fold",                   1, "number of encodings per time step batch")
    .item("restart",               "", "a checkpoint file to restart from")
    .item("outdir",               ".", "the directory where output files are written");
}

inline auto make_solver_data(const mara::config_t& cfg)
{
    return solver_data_t{
        cfg.get_int("block_size"),
        cfg.get_int("depth"),
        cfg.get_double("omega_frame"),
        cfg.get_double("mach_number"),
        cfg.get_double("nu"),
        cfg.get_double("domain_radius"),
        cfg.get_double("softening_length"),
        cfg.get_double("eccentricity"),
        cfg.get_double("mass_ratio"),
        cfg.get_double("sink_radius"),
        cfg.get_double("sink_rate"),
        cfg.get_double("buffer_rate"),
        cfg.get_double("buffer_scale"),
        cfg.get_double("cfl"),
    };
}




//=============================================================================
unit_specific_energy sound_speed_squared(numeric::array_t<unit_length, 2> p, unit_length softening_length, unit_scalar mach_number);

nd::shared_array<numeric::array_t<unit_length, 2>, 2> face_coordinates(int block_size, unit_length domain_radius, bsp::tree_index_t<2> block, bsp::uint axis);
nd::shared_array<numeric::array_t<unit_length, 2>, 2> cell_coordinates(int block_size, unit_length domain_radius, bsp::tree_index_t<2> block);

nd::shared_array<unit_rate, 2>                             buffer_rate_field_array(solver_data_t solver_data, bsp::tree_index_t<2> block);
nd::shared_array<numeric::array_t<unit_acceleration, 2>, 2> centrifugal_term_array(solver_data_t solver_data, bsp::tree_index_t<2> block);
nd::shared_array<iso2d::conserved_density_t, 2>            initial_conserved_array(solver_data_t solver_data, bsp::tree_index_t<2> block);




//=============================================================================
unit_length       smallest_cell_size(solver_data_t solver_data);
primitive_array_t recover_primitive_array(conserved_array_t uc);
primitive_array_t estimate_gradient(primitive_array_t pc, bsp::uint axis, double theta);
conserved_array_t updated_conserved(
    conserved_array_t uc,
    primitive_array_t pe,
    unit_time time,
    unit_time dt,
    bsp::tree_index_t<2> block,
    solver_data_t solver_data);

} // namespace minidisk
