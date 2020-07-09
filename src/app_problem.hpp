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
#include "core_rational.hpp"
#include "parallel_computable.hpp"




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
     * @brief      Run a set of side effects on the solution, using the
     *             recurrence rules set in the run config, and the last
     *             invocation time in the given schedule. Update the schedule
     *             with the new last invocation time.
     *
     * @param[in]  run_config  The run config
     * @param[in]  schedule    The schedule instance
     * @param[in]  solution    The solution instance
     */
    void side_effects(const mara::config_t& run_config, schedule_t& schedule, std::any solution, mpr::execution_monitor_t monitor) const;




    //=============================================================================
    void print_iteration_message(const mara::config_t& run_config, const schedule_t& schedule, std::any solution, mpr::execution_monitor_t monitor) const;
    void write_checkpoint(const mara::config_t& run_config, schedule_t& schedule, std::any solution) const;
    void print_task_graph(const mara::config_t& run_config, schedule_t& schedule, std::any solution) const;



    /**
     * @brief      Generate initial data by examining the run config for a
     *             restart file, invoking the read_solution virtual function if
     *             one is given, but otherwise returning the given initial
     *             solution data.
     *
     * @param[in]  run_config               The run config
     * @param[in]  solution_if_not_restart  The solution to return if no restart
     *                                      file was given
     *
     * @return     A type-erased solution struct
     */
    std::any initial_solution(const mara::config_t& run_config, std::any solution_if_not_restart) const;




    /**
     * @brief      Get the time from solution. Your derived class must implement
     *             this method to return the simulation time from the solution
     *             struct.
     *
     * @param[in]  solution  The (type-erased) solution instance
     *
     * @return     The time
     */
    virtual dimensional::unit_time get_time_from_solution(std::any solution) const = 0;
    virtual rational::number_t get_iteration_from_solution(std::any solution) const = 0;
    virtual unsigned long get_zone_count(std::any solution) const = 0;




    /**
     * @brief      Write a solution to a checkpoint file. Your derived class
     *             must implement this method to write the given solution struct
     *             to a checkpoint file with the specified filename. This method
     *             is invoked sequentially on each MPI process, so you should
     *             not do any MPI barriers in your implementation; just open the
     *             file, write to it, and close it.
     *
     * @param[in]  filename  The checkpoint filename
     * @param[in]  solution  The (type-erased) solution struct
     */
    virtual void write_solution(std::string filename, std::any solution) const = 0;




    /**
     * @brief      Read a solution from a checkpoint file. Your derived class
     *             must implement this method to read a solution struct from a
     *             checkpoint with the specified filename.
     *
     * @param[in]  filename  The checkpoint filename
     *
     * @return     A type-erased solution struct
     */
    virtual std::any read_solution(std::string filename) const = 0;




    /**
     * @brief      Return a set of computable nodes in the solution. This set is
     *             used to print the task graph in graphviz format, and in
     *             driving the evaluation of computable graphs.
     *
     * @param[in]  solution  The (type-erased) solution struct
     *
     * @return     A set of computable nodes
     */
    virtual mpr::node_set_t get_computable_nodes(std::any solution) const = 0;




    /**
     * @brief      This method must be over-ridden to return a module name. The
     *             module name is used to infer the problem type from restart
     *             files, and to interpret data in plotting scripts. It is
     *             probably synonymous with the type of solution struct being
     *             used. The string is written to HDF5 checkpoint files as
     *             'module' in the root group.
     *
     * @return     A string identifying the module name.
     */
    virtual std::string get_module_name() const = 0;
};




/**
 * @brief      Convenience method to invoke the problem initial_solution member
 *             function, handling the type-erasure automatically.
 *
 * @param[in]  problem                  The problem instance
 * @param[in]  run_config               The run config
 * @param[in]  solution_if_not_restart  The solution to return if no restart
 *                                      file was given
 *
 * @tparam     SolutionType             The type of the solution struct
 *
 * @return     A solution struct
 */
template<typename SolutionType>
SolutionType initial_solution(const problem_base_t& problem, const mara::config_t& run_config, SolutionType solution_if_not_restart)
{
    return std::any_cast<SolutionType>(problem.initial_solution(run_config, std::move(solution_if_not_restart)));
}




/**
 * @brief      Generate an initial schedule instance, either by reading it from
 *             a checkpoint file if a restart file is given in the run_config,
 *             or otherwise the one given (defaulting to a fresh schedule).
 *
 * @param[in]  run_config               The run config
 * @param[in]  schedule_if_not_restart  The schedule to return if no restart
 *                                      file was given
 *
 * @return     A schedule struct
 */
schedule_t initial_schedule(const mara::config_t& run_config, schedule_t schedule_if_not_restart={});




/**
 * @brief      Read and return the run_config group of a restart file if one is
 *             given in the arguments as e.g. restart=chkpt.1234.h5, otherwise
 *             return an empty map.
 *
 * @param[in]  args  The command line arguments
 *
 * @return     A configuration parameter map
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
