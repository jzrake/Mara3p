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




#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_rational.hpp"
#include "core_bqo_tree.hpp"
#include "core_util.hpp"
#include "minidisk_io.hpp"




//=============================================================================
void h5::read(const Group& group, std::string name, minidisk::side_effect_t& side_effect)
{
    auto sgroup = group.open_group(name);
    read(sgroup, "next_due", side_effect.next_due);
    read(sgroup, "count", side_effect.count);
}

void h5::read(const Group& group, std::string name, minidisk::schedule_t& schedule)
{
    auto sgroup = group.open_group(name);
    read(sgroup, "time_series", schedule.time_series);
    read(sgroup, "checkpoint", schedule.checkpoint);
}

void h5::read(const Group& group, std::string name, minidisk::solution_t& solution)
{
    auto sgroup = group.open_group(name);
    auto ugroup = sgroup.open_group("conserved");

    read(sgroup, "iteration", solution.iteration);
    read(sgroup, "time", solution.time);

    for (auto block_name : ugroup)
    {
        auto block = bsp::read_tree_index<2>(block_name);
        auto value = read<minidisk::conserved_array_t>(ugroup, block_name);
        solution.conserved = insert(solution.conserved, block, value);
    }
}




//=============================================================================
void h5::write(const Group& group, std::string name, const minidisk::side_effect_t& side_effect)
{
    auto sgroup = group.require_group(name);
    write(sgroup, "next_due", side_effect.next_due);
    write(sgroup, "count", side_effect.count);
}

void h5::write(const Group& group, std::string name, const minidisk::schedule_t& schedule)
{
    auto sgroup = group.require_group(name);
    write(sgroup, "time_series", schedule.time_series);
    write(sgroup, "checkpoint", schedule.checkpoint);
}

void h5::write(const Group& group, std::string name, const minidisk::solution_t& solution)
{
    auto sgroup = group.require_group(name);
    auto ugroup = sgroup.require_group("conserved");

    write(sgroup, "iteration", solution.iteration);
    write(sgroup, "time", solution.time);

    sink(indexify(solution.conserved), util::apply_to([&ugroup] (auto index, auto value)
    {
        write(ugroup, bsp::format_tree_index(index), value);
    }));
}
