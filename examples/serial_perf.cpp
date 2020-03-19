#include <cstdio>
#include <chrono>
#include "core_numeric_array.hpp"
#include "core_numeric_tuple.hpp"
#include "app_serial_std_vector.hpp"
#include "app_serial_ndarray.hpp"
#include "app_serial_numeric_tuple.hpp"

#define R 1000
#define N 256




//=============================================================================
template<typename T>
auto time_serialize(const T& value)
{
    auto start = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < R; ++i)
    {
        serial::dumps(value);
    }
    return std::chrono::high_resolution_clock::now() - start;
}




//=============================================================================
int main()
{
    auto P0 = std::vector<double>(N * N * 3);
    auto P1 = nd::zeros(N, N) | nd::map([] (auto) { return numeric::array(1.0, 2.0, 3.0); }) | nd::to_shared();
    auto P2 = nd::zeros(N, N) | nd::map([] (auto) { return numeric::tuple(1.0, 2.0, 3.0); }) | nd::to_shared();

    std::printf("std::vector ........ %lfs (%ld bytes)\n", 1e-9 * time_serialize(P0).count(), serial::dumps(P0).size());
    std::printf("numeric::array ..... %lfs (%ld bytes)\n", 1e-9 * time_serialize(P1).count(), serial::dumps(P1).size());
    std::printf("numeric::tuple ..... %lfs (%ld bytes)\n", 1e-9 * time_serialize(P2).count(), serial::dumps(P2).size());

    return 0;
}
