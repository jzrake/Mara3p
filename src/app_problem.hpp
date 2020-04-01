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
#include <any>
#include <string>
#include "app_config.hpp"
#include "core_dimensional.hpp"




namespace mara
{




/**
 * @brief      Data structure holding recurrence information for a single
 *             side-effect.
 */
struct side_effect_t
{
    dimensional::unit_time next_due = 0.0;
    int count = 0;
};




/**
 * @brief      Data structure holding various types of side effects a problem
 *             might issue.
 */
struct schedule_t
{
    side_effect_t checkpoint;
    side_effect_t time_series;
};




/**
 * @brief      This class describes a problem setup.
 */
class problem_base_t
{
public:




    /**
     * @brief      { function_description }
     *
     * @param[in]  run_config  The run configuration
     * @param[in]  time        The time
     * @param[in]  schedule    The schedule
     * @param[in]  solution    The solution
     *
     * @return     { description_of_the_return_value }
     */
    schedule_t side_effects(const mara::config_t& run_config, schedule_t schedule, std::any solution) const;




    /**
     * @brief      Gets the time from solution.
     *
     * @param[in]  solution  The solution
     *
     * @return     The time from solution.
     */
    virtual dimensional::unit_time get_time_from_solution(std::any solution) const = 0;




    /**
     * @brief      Writes a solution to checkpoint.
     *
     * @param[in]  filename  The filename
     * @param[in]  solution  The solution
     */
    virtual void write_solution(std::string filename, std::any solution) const = 0;




    /**
     * @brief      Reads a solution from checkpoint.
     *
     * @param[in]  filename  The filename
     *
     * @return     { description_of_the_return_value }
     */
    virtual std::any read_solution(std::string filename) const = 0;
};




/**
 * @brief      { function_description }
 *
 * @param[in]  cfg                      The configuration
 * @param[in]  solution_if_not_restart  The solution if not restart
 *
 * @return     { description_of_the_return_value }
 */
std::any initial_solution(const mara::config_t& cfg, std::any solution_if_not_restart);




/**
 * @brief      { function_description }
 *
 * @param[in]  cfg                      The configuration
 * @param[in]  schedule_if_not_restart  The schedule if not restart
 *
 * @return     { description_of_the_return_value }
 */
schedule_t initial_schedule(const mara::config_t& cfg, schedule_t schedule_if_not_restart={});




/**
 * @brief      { function_description }
 *
 * @param[in]  args  The arguments
 *
 * @return     { description_of_the_return_value }
 */
config_parameter_map_t restart_run_config(const config_string_map_t& args);

} // namespace mara




//=============================================================================
namespace h5 {

class Group;

void read(const Group& group, std::string name, mara::side_effect_t& side_effect);
void read(const Group& group, std::string name, mara::schedule_t& schedule);

void write(const Group& group, std::string name, const mara::side_effect_t& side_effect);
void write(const Group& group, std::string name, const mara::schedule_t& schedule);

} // namespace h5
