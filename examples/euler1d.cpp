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




//=============================================================================
#include <chrono>
#include <tuple>

template<typename Function, typename... Args>
auto time_execution(Function&& func, Args&&... args)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto result = std::forward<Function>(func)(std::forward<Args>(args)...);
    auto stop = std::chrono::high_resolution_clock::now();
    auto ms = 1e-6 * std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    return std::make_pair(std::move(result), ms);
};

template<typename T>
auto construct()
{
    return [] (auto x) { return T(x); };
}




//=============================================================================
#include "app_hdf5.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "physics_euler.hpp"




//=============================================================================
static const auto gamma_law_index = 5. / 3;
static const auto xhat = geometric::unit_vector_t{{1.0, 0.0, 0.0}};
static const auto num_cells = 8192;




//=============================================================================
namespace nd {

inline auto to_shared_prof(const char* message)
{
    return [/*message*/] (auto x)
    {
        return to_shared(x);
        // auto result = time_execution([] (auto x) { return nd::to_shared(x); }, x);
        // std::printf("\t[%-20s]: %5.4e\n", message, result.second);
        // return result.first;
    };
}

}
    




//=============================================================================
namespace mesh {

inline auto vertices(nd::uint cell_count)
{
    return nd::linspace(0.0, 1.0, cell_count + 1) | nd::map(construct<dimensional::unit_length>());
}

inline auto cell_centers(nd::uint cell_count)
{
    return vertices(cell_count) | nd::adjacent_mean();
}

inline auto cell_spacings(nd::uint cell_count)
{
    return vertices(cell_count) | nd::adjacent_diff();
}

} // namespace mesh




//=============================================================================
euler::primitive_t initial_condition(dimensional::unit_length x)
{
    return x < dimensional::unit_length(0.5) ? euler::primitive(1.0, 1.0) : euler::primitive(0.125, 0.1);
}

auto riemann_solver_for(geometric::unit_vector_t nhat)
{
    return [nhat] (auto t)
    {
        return euler::riemann_hlle(std::get<0>(t), std::get<1>(t), nhat, gamma_law_index);
    };
}




//=============================================================================
struct state_t
{
    dimensional::unit_time time = 0.0;
    rational::number_t iteration = 0;
    nd::shared_array<euler::conserved_density_t, 1> conserved;
};

state_t initial_state()
{
    auto p0 = mesh::cell_centers(num_cells) | nd::map(initial_condition);
    auto u0 = p0 | nd::map([] (auto p) { return euler::conserved_density(p, gamma_law_index); });

    return {0.0, 0, u0 | nd::to_shared()};
}

state_t advance(state_t state)
{
    auto dt = dimensional::unit_time(0.1 / num_cells);
    auto dx = mesh::cell_spacings(num_cells);
    auto u0 = state.conserved;
    auto p0 = u0 | nd::map([] (auto u) { return euler::recover_primitive(u, gamma_law_index); });
    auto f0 = p0 | nd::extend_zero_gradient() | nd::adjacent_zip() | nd::map(riemann_solver_for(xhat)) | nd::to_shared_prof("Riemann");
    auto df = f0 | nd::adjacent_diff();
    auto u1 = u0 - df * dt / dx | nd::to_shared_prof("U + dU");

    return {
        state.time + dt,
        state.iteration + rational::number(1),
        u1,
    };
}




//=============================================================================
int main()
{
    auto state = initial_state();
    auto total_kzps = 0.0;

    while (state.time < dimensional::unit_time(0.3))
    {
        auto [next_state, ms] = time_execution(advance, state);
        std::printf("[%04lu] t=%.4e kzps=%3.2lf\n", long(state.iteration), state.time.value, double(num_cells) / ms);

        total_kzps += double(num_cells) / ms;
        state = next_state;
    }

    std::printf("finished: mean kzps = %3.2lf\n", total_kzps / double(state.iteration));

    auto file = h5::File("euler.h5", "w");

    h5::write(file, "primitive",
        state.conserved
        | nd::map([] (auto u) { return euler::recover_primitive(u, gamma_law_index); })
        | nd::to_shared());

    h5::write(file, "cell_centers", mesh::cell_centers(num_cells) | nd::map([] (auto x) { return x.value; }) | nd::to_shared());

    return 0;
}
