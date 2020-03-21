#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include "app_config.hpp"
#include "parallel_mpi.hpp"
#include "parallel_thread_pool.hpp"




//=============================================================================
double compute_pi(int seed, int num_samples)
{
    auto generator = std::default_random_engine(seed);
    auto distribution = std::uniform_real_distribution<double>(-1.0, 1.0);
    auto total = 0.0;

    for (int i = 0; i < num_samples; ++i)
    {
        auto x = distribution(generator);
        auto y = distribution(generator);
        total += x * x + y * y < 1.0;
    }
    return 4.0 * total / num_samples;
}




//=============================================================================
int main(int argc, const char* argv[])
{
    auto mpi_session  = mpi::Session();
    auto cfg_template = mara::config_template()
    .item("threads",              1)   // number of threads to use (zero for hardware concurrency)
    .item("num_samples",        1e5)   // number of Monte-Carlo samples per batch
    .item("num_batches",        1e3);  // number of total batches (not per MPI process)

    auto cfg = cfg_template.create().update(mara::argv_to_string_map(argc, argv));
    auto num_threads = cfg.get_int("threads");
    auto thread_pool = mara::ThreadPool(num_threads == 0 ? std::thread::hardware_concurrency() : num_threads);
    auto num_samples = int(cfg.get_double("num_samples"));
    auto num_batches = int(cfg.get_double("num_batches") / mpi::comm_world().size());
    auto results     = std::vector<std::future<double>>();




    //=========================================================================
    if (mpi::comm_world().rank() == 0)
    {
        mara::pretty_print(std::cout, "config", cfg);
    }

    for (int i = 0; i < num_batches; ++i)
    {
        auto seed = i + 100000 * mpi::comm_world().rank();
        results.push_back(thread_pool.enqueue(compute_pi, seed, num_samples));
    }




    //=========================================================================
    auto start = std::chrono::high_resolution_clock::now();
    auto total = 0.0;

    for (auto& result : results)
    {
        total += result.get();
    }
    auto all_total = mpi::comm_world().reduce(0, total, mpi::operation::sum);
    auto finish = std::chrono::high_resolution_clock::now();




    //=========================================================================
    if (mpi::comm_world().rank() == 0)
    {
        std::printf("MPI processes .................. %d\n", mpi::comm_world().size());
        std::printf("threads / MPI process .......... %lu\n", thread_pool.size());
        std::printf("total compute units ............ %lu\n", mpi::comm_world().size() * thread_pool.size());
        std::printf("compute-seconds ................ %lf\n", mpi::comm_world().size() * thread_pool.size() * (finish - start).count() * 1e-9);
        std::printf("pi estimate .................... %lf\n", all_total / num_batches / mpi::comm_world().size());
    }
    return 0;
}
