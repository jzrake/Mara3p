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




#pragma once
#include "app_hdf5.hpp"
#include "core_ndarray.hpp"
#include "core_numeric_tuple.hpp"
#include "core_numeric_array.hpp"
#include "core_dimensional.hpp"




//=============================================================================
namespace h5 {




//=============================================================================
template<nd::uint Rank, long... N1, long... N2, long... N3, unsigned long... D1, unsigned long... D2, unsigned long... D3>
struct h5::hdf5_conversion_to_hdf5_writable<nd::shared_array<numeric::tuple_t<dimensional::quantity_t<N1, N2, N3, D1, D2, D3>...>, Rank>>
{
    using type = nd::shared_array<numeric::array_t<double, sizeof...(N1)>, Rank>;

    auto operator()(const nd::shared_array<numeric::tuple_t<dimensional::quantity_t<N1, N2, N3, D1, D2, D3>...>, 1>& value) const
    {
        return map(value, [] (auto q)
        {
            return std::apply([] (auto... t)
            {
                return numeric::array(t.value...); }, q.impl);
            }) | nd::to_shared();
    }
};

template<nd::uint Rank, long... N1, long... N2, long... N3, unsigned long... D1, unsigned long... D2, unsigned long... D3>
struct h5::hdf5_conversion_from_hdf5_writable<nd::shared_array<numeric::tuple_t<dimensional::quantity_t<N1, N2, N3, D1, D2, D3>...>, Rank>>
{
    using type = nd::shared_array<numeric::array_t<double, sizeof...(N1)>, Rank>;

    auto operator()(const type& value) const
    {
        return map(value, [] (auto x)
        {
            return std::apply([] (auto... xs)
            {
                return numeric::tuple(dimensional::quantity_t<N1, N2, N3, D1, D2, D3>{xs}...);
            }, x.impl);
        }) | nd::to_shared();
    }
};

} // namespace h5




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "app_hdf5_ndarray.hpp"
#include "core_unit_test.hpp"




//=============================================================================
inline void test_hdf5_ndarray_dimensional()
{
    auto test_read_write = [] (auto value)
    {
        {
            auto file = h5::File("test.h5", h5::File::Access::Truncate);
            h5::write(file, "value", value);
        }
        {
            auto file = h5::File("test.h5", h5::File::Access::Read);
            require(all(h5::read<decltype(value)>(file, "value") == value));
        }
    };

    test_read_write(nd::from(
        numeric::tuple(1.2 * dimensional::unit_mass{}, 1.2 * dimensional::unit_energy{}),
        numeric::tuple(1.3 * dimensional::unit_mass{}, 1.3 * dimensional::unit_energy{})) | nd::to_shared());
}

#endif // DO_UNIT_TESTS
