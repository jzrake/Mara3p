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
#include "core_ndarray.hpp"
#include "physics_mhd.hpp"




//=============================================================================
namespace mara {




//=============================================================================
struct mhd_boundary_extension
{
    using primitive_array = nd::shared_array<mhd::primitive_t, 3>;
    using magnetic_array  = nd::shared_array<mhd::unit_magnetic_field, 3>;

    virtual ~mhd_boundary_extension() {}
    virtual primitive_array extend_primitive(primitive_array pc) const = 0;
    virtual magnetic_array extend_magnetic1(magnetic_array bf1) const = 0;
    virtual magnetic_array extend_magnetic2(magnetic_array bf2) const = 0;
    virtual magnetic_array extend_magnetic3(magnetic_array bf3) const = 0;
};

std::shared_ptr<mhd_boundary_extension> local_periodic_boundary_extension();




//=============================================================================
nd::shared_array<mhd::primitive_t, 3> primitive_array(
    nd::shared_array<mhd::conserved_density_t, 3> uc,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf1,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf2,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf3,
    double gamma_law_index);




//=============================================================================
std::tuple<
    nd::shared_array<mhd::flux_vector_t, 3>,
    nd::shared_array<mhd::flux_vector_t, 3>,
    nd::shared_array<mhd::flux_vector_t, 3>,
    nd::shared_array<mhd::unit_electric_field, 3>,
    nd::shared_array<mhd::unit_electric_field, 3>,
    nd::shared_array<mhd::unit_electric_field, 3>>
flux_arrays(
    nd::shared_array<mhd::primitive_t, 3> pc,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf1,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf2,
    nd::shared_array<mhd::unit_magnetic_field, 3> bf3);




//=============================================================================
std::tuple<
    dimensional::unit_time,
    nd::shared_array<mhd::conserved_density_t, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>>
advance(dimensional::unit_time time,
        nd::shared_array<mhd::conserved_density_t, 3> uc,
        nd::shared_array<mhd::unit_magnetic_field, 3> bf1,
        nd::shared_array<mhd::unit_magnetic_field, 3> bf2,
        nd::shared_array<mhd::unit_magnetic_field, 3> bf3,
        dimensional::unit_length dl,
        const mhd_boundary_extension& boundary_extension);




//=============================================================================
std::tuple<
    nd::shared_array<mhd::conserved_density_t, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>,
    nd::shared_array<mhd::unit_magnetic_field, 3>>
construct_conserved(
    std::function<mhd::primitive_t(geometric::euclidean_vector_t<dimensional::unit_length>, mhd::magnetic_field_vector_t)> primitive,
    std::function<mhd::vector_potential_t(geometric::euclidean_vector_t<dimensional::unit_length>)> vector_potential,
    nd::uint block_size);

} // namespace mara
