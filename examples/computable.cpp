#include <iostream>
#include "parallel_mpi.hpp"
#include "app_serial_std_tuple.hpp"
#include "app_serial_std_string.hpp"
#include "parallel_computable.hpp"




//=============================================================================
auto concat(mpr::computable<std::string> a, mpr::computable<std::string> b)
{
    return zip(a, b).name("zip") | mpr::mapv([] (auto a, auto b) { return a + b; });
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
    auto abbb = concat(ab, bb).name("abbb");
    auto baaa = concat(ba, aa).name("baaa");
    auto baab = concat(ba, ab).name("baab");
    auto baba = concat(ba, ba).name("baba");
    auto babb = concat(ba, bb).name("babb");
    auto bbaa = concat(bb, aa).name("bbaa");
    auto bbab = concat(bb, ab).name("bbab");
    auto bbbb = concat(bb, bb).name("bbbb");


    auto master = concat(aaaa, aaab, aaba, aabb, abaa, abab, abbb, baaa, baab, baba, babb, bbaa, bbab, bbbb).name("master");
    mpr::compute_mpi(master);

    if (master.has_value())
    {
        std::cout << "master.value = " << master.value() << std::endl;
    }

    MPI_Finalize();
    return 0;
}
