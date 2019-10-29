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




#include "app_hdf5.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "core_ndarray.hpp"
#include "physics_euler.hpp"




//=============================================================================
static const auto gamma_law_index = 5. / 3;
static const auto xhat = geometric::unit_vector_t{{1.0, 0.0, 0.0}};
static const auto num_cells = 10000;




//=============================================================================
namespace nd {

inline auto adjacent_zip(uint axis=0)
{
    return [axis] (auto x)
    {
        return zip(select(x, axis, 0, -1), select(x, axis, 1));
    };
}

inline auto adjacent_mean(uint axis=0)
{
    return [axis] (auto x)
    {
        return 0.5 * (select(x, axis, 0, -1) + select(x, axis, 1));
    };
}

inline auto adjacent_diff(uint axis=0)
{
    return [axis] (auto x)
    {
        return select(x, axis, 1) - select(x, axis, 0, -1);
    };
}

inline auto extend_periodic(uint axis=0)
{
    return [axis] (auto x)
    {
        return select(x, axis, -1) | nd::concat(x, axis) | nd::concat(select(x, axis, 0, 1), axis);
    };
}

}




//=============================================================================
template<typename T>
auto construct()
{
    return [] (auto x) { return T(x); };
}



#include <chrono>

template<typename Function, typename... Args>
auto time_execution(Function&& func, Args&&... args)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto result = std::forward<Function>(func)(std::forward<Args>(args)...);
    auto stop = std::chrono::high_resolution_clock::now();
    return std::make_pair(std::move(result), stop - start);
};




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
    auto f0 = p0 | nd::extend_periodic() | nd::adjacent_zip() | nd::map(riemann_solver_for(xhat));
    auto df = f0 | nd::adjacent_diff();
    auto u1 = u0 - df * dt / dx;

    return {
        state.time + dt,
        state.iteration + rational::number(1),
        u1 | nd::to_shared()
    };
}




//=============================================================================
int main()
{
    auto state = initial_state();


    while (state.iteration < 1000)
    {
        auto [next_state, duration] = time_execution(advance, state);

        auto ms = 1e-6 * std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
        auto kzps = double(num_cells) / ms;

        state = next_state;
        std::printf("[%04lu] t=%4.4lf kzps=%3.2lf\n", state.iteration.num, state.time.value, kzps);
    }


    auto file = h5::File("euler.h5", "w");

    h5::write(file, "primitive",
        state.conserved
        | nd::map([] (auto u) { return euler::recover_primitive(u, gamma_law_index); })
        | nd::to_shared());

    h5::write(file, "cell_centers", mesh::cell_centers(num_cells) | nd::map([] (auto x) { return x.value; }) | nd::to_shared());

    return 0;
}
