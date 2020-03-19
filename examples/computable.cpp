#include <iostream>
#include <string>
#include "app_serial.hpp"
#include "app_serial_std_vector.hpp"
#include "parallel_computable.hpp"
#include "parallel_mpi.hpp"

#define DUMP_DELAY std::chrono::milliseconds(50)
#define OPERATION_DELAY std::chrono::milliseconds(100)




//=============================================================================
template<>
struct serial::conversion_to_serializable_t<std::string>
{
    using type = std::vector<char>;
    auto operator()(std::string value) const
    {
        std::this_thread::sleep_for(DUMP_DELAY);
        return std::vector<char>(value.begin(), value.end());
    }
};

template<>
struct serial::conversion_from_serializable_t<std::string>
{
    using type = std::vector<char>;
    auto operator()(std::vector<char> value) const { return std::string(value.begin(), value.end()); }
};

template<>
struct serial::is_serializable_t<std::string> : std::true_type {};




//=============================================================================
auto concat(mpr::computable<std::string> a, mpr::computable<std::string> b)
{
    return zip(a, b).name("zip") | mpr::mapv([] (auto a, auto b)
    {
        std::this_thread::sleep_for(OPERATION_DELAY);
        return a + b;
    });
}

template<typename T, typename U, typename... Vs>
auto concat(mpr::computable<T> t, mpr::computable<U> u, mpr::computable<Vs>... vs)
{
    return concat(concat(t, u).name("concat"), vs...);
}




//=============================================================================
int main()
{
    MPI_Init(nullptr, nullptr);

    auto a = mpr::just<std::string>("a").name("a");
    auto b = mpr::just<std::string>("b").name("b");
    auto aa = concat(a, a).name("aa");
    auto ab = concat(a, b).name("ab");
    auto ba = concat(b, a).name("ba");
    auto bb = concat(b, b).name("bb");

    auto aaaa = concat(aa, aa).name("aaaa");
    auto aaab = concat(aa, ab).name("aaab");
    auto aaba = concat(aa, ba).name("aaba");
    auto aabb = concat(aa, bb).name("aabb");
    auto abaa = concat(ab, aa).name("abaa");
    auto abab = concat(ab, ab).name("abab");
    auto abba = concat(ab, ba).name("abaa");
    auto abbb = concat(ab, bb).name("abbb");
    auto baaa = concat(ba, aa).name("baaa");
    auto baab = concat(ba, ab).name("baab");
    auto baba = concat(ba, ba).name("baba");
    auto babb = concat(ba, bb).name("babb");
    auto bbaa = concat(bb, aa).name("bbaa");
    auto bbab = concat(bb, ab).name("bbab");
    auto bbba = concat(bb, ba).name("bbba");
    auto bbbb = concat(bb, bb).name("bbbb");

    mpr::compute_mpi(aaaa, aaab, aaba, aabb, abaa, abab, abba, abbb, baaa, baab, baba, babb, bbaa, bbab, bbba, bbbb);

    mpi::comm_world().invoke([&] ()
    {
        if (aaaa.has_value()) { printf("%s compued on rank %d: %s\n", "aaaa", mpi::comm_world().rank(), aaaa.value().data()); }
        if (aaab.has_value()) { printf("%s compued on rank %d: %s\n", "aaab", mpi::comm_world().rank(), aaab.value().data()); }
        if (aaba.has_value()) { printf("%s compued on rank %d: %s\n", "aaba", mpi::comm_world().rank(), aaba.value().data()); }
        if (aabb.has_value()) { printf("%s compued on rank %d: %s\n", "aabb", mpi::comm_world().rank(), aabb.value().data()); }
        if (abaa.has_value()) { printf("%s compued on rank %d: %s\n", "abaa", mpi::comm_world().rank(), abaa.value().data()); }
        if (abab.has_value()) { printf("%s compued on rank %d: %s\n", "abab", mpi::comm_world().rank(), abab.value().data()); }
        if (abba.has_value()) { printf("%s compued on rank %d: %s\n", "abba", mpi::comm_world().rank(), abba.value().data()); }
        if (abbb.has_value()) { printf("%s compued on rank %d: %s\n", "abbb", mpi::comm_world().rank(), abbb.value().data()); }
        if (baaa.has_value()) { printf("%s compued on rank %d: %s\n", "baaa", mpi::comm_world().rank(), baaa.value().data()); }
        if (baab.has_value()) { printf("%s compued on rank %d: %s\n", "baab", mpi::comm_world().rank(), baab.value().data()); }
        if (baba.has_value()) { printf("%s compued on rank %d: %s\n", "baba", mpi::comm_world().rank(), baba.value().data()); }
        if (babb.has_value()) { printf("%s compued on rank %d: %s\n", "babb", mpi::comm_world().rank(), babb.value().data()); }
        if (bbaa.has_value()) { printf("%s compued on rank %d: %s\n", "bbaa", mpi::comm_world().rank(), bbaa.value().data()); }
        if (bbab.has_value()) { printf("%s compued on rank %d: %s\n", "bbab", mpi::comm_world().rank(), bbab.value().data()); }
        if (bbba.has_value()) { printf("%s compued on rank %d: %s\n", "bbba", mpi::comm_world().rank(), bbba.value().data()); }
        if (bbbb.has_value()) { printf("%s compued on rank %d: %s\n", "bbbb", mpi::comm_world().rank(), bbbb.value().data()); }
    });

    MPI_Finalize();
    return 0;
}
