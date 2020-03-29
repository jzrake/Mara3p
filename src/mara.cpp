/**
 ==============================================================================
 Copyright 2019 - 2020, Jonathan Zrake

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




#include "core_unit_test.hpp"
#include "module_euler.hpp"
#include "app_control.hpp"
#include "app_hdf5.hpp"

using namespace std::placeholders;




//=============================================================================
int test_app();
int test_core();
int test_mesh();
int test_model();
int test_physics();




//=============================================================================
static void run_all_tests()
{
    start_unit_tests();
    test_app();
    test_core();
    test_mesh();
    test_model();
    test_physics();
    report_test_results();
}




//=============================================================================
static void run_euler2d()
{
    using namespace modules::euler2d;

    auto shocktube_2d = [] (auto x)
    {
        if (sqrt(sum(x * x)) < dimensional::unit_length(0.2))
        {
            return euler::primitive(1.0, 0.0, 0.0, 0.0, 1.0);
        }
        else
        {
            return euler::primitive(0.1, 0.0, 0.0, 0.0, 0.125);
        }
    };

    auto tfinal           = dimensional::unit_time(0.0);
    auto plm_theta        = 2.0;
    auto gamma_law_index  = 5. / 3;
    auto geom             = mesh_geometry_t();
    auto mesh             = bsp::quadtree([geom] (auto block)
    {
        auto p = geom.block_centroid(block);
        auto r = sqrt(sum(p * p));
        return r < unit_length(1.0) / double(block.level);
    }, 4);

    auto u                = mpr::compute_all(initial_conserved_tree(mesh, geom, shocktube_2d, gamma_law_index));
    auto solution         = solution_t{0, 0.0, u};
    auto dx = smallest_cell_size(mesh, geom);
    auto kz = total_cells(mesh, geom) / 1e3;
    auto vm = dimensional::unit_velocity(10.0);


    while (solution.time < tfinal)
    {
        auto dt = dx / vm;
        auto scheme = std::bind(updated_solution, _1, dt, geom, plm_theta, gamma_law_index);

        solution = control::advance_runge_kutta(scheme, 2, solution);

        auto start = std::chrono::high_resolution_clock::now();
        mpr::compute_all(solution.conserved, 0);
        auto ticks = std::chrono::high_resolution_clock::now() - start;
        std::printf("[%06ld] t=%.4lf kzps=%.2lf\n", long(solution.iteration), solution.time.value, kz / (ticks.count() * 1e-9));
    }

    auto h5f = h5::File("euler2d.h5", "w");

    h5::write(h5f, "solution", solution);
    h5::write(h5f, "module", std::string("euler2d"));
}




//=============================================================================
int main(int argc, const char* argv[])
{
    if (argc == 1)
    {
        std::printf("usage: mara [command] <key-val parameters>\n");
        return 0;
    }

    auto command = std::string(argv[1]);

    if (false) {}
    else if (command == "test")    run_all_tests();
    else if (command == "euler2d") run_euler2d();
    else std::printf("unknown command: %s\n", command.data());

    return 0;
}
