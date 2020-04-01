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




#include "app_problem.hpp"
#include "app_hdf5_config.hpp"
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_rational.hpp"
#include "core_util.hpp"




//=============================================================================
mara::schedule_t mara::initial_schedule(const mara::config_t& cfg, schedule_t schedule_if_not_restart)
{
    auto restart = cfg.get_string("restart");

    if (! restart.empty())
    {
        auto result = schedule_t{};
        auto h5f = h5::File(restart, "r");
        read(h5f, "schedule", result);
        return result;
    }
    return schedule_t();
}

mara::config_parameter_map_t mara::restart_run_config(const mara::config_string_map_t& args)
{
    if (args.count("restart"))
    {
        auto file = h5::File(args.at("restart"), "r");
        return h5::read<mara::config_parameter_map_t>(file, "run_config");
    }
    return mara::config_parameter_map_t{};
}




//=============================================================================
void h5::read(const Group& group, std::string name, mara::side_effect_t& side_effect)
{
    auto sgroup = group.open_group(name);
    read(sgroup, "next_due", side_effect.next_due);
    read(sgroup, "count", side_effect.count);
}

void h5::read(const Group& group, std::string name, mara::schedule_t& schedule)
{
    auto sgroup = group.open_group(name);
    read(sgroup, "time_series", schedule.time_series);
    read(sgroup, "checkpoint", schedule.checkpoint);
}




//=============================================================================
void h5::write(const Group& group, std::string name, const mara::side_effect_t& side_effect)
{
    auto sgroup = group.require_group(name);
    write(sgroup, "next_due", side_effect.next_due);
    write(sgroup, "count", side_effect.count);
}

void h5::write(const Group& group, std::string name, const mara::schedule_t& schedule)
{
    auto sgroup = group.require_group(name);
    write(sgroup, "time_series", schedule.time_series);
    write(sgroup, "checkpoint", schedule.checkpoint);
}
