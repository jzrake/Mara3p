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
#include "core_ndarray.hpp"
#include "app_hdf5.hpp"




//=============================================================================
namespace h5 {




//=============================================================================
template<typename ProviderType, nd::uint Rank>
struct hdf5_container_address<nd::array_t<ProviderType, Rank>>
{
    const void* operator()(const nd::array_t<ProviderType, Rank>& value)
    {
        return value.data();
    }
    void* operator()(nd::array_t<ProviderType, Rank>& value)
    {
        return value.data();
    }
};

template<typename ProviderType, nd::uint Rank>
struct hdf5_dataspace_creation<nd::array_t<ProviderType, Rank>>
{
    auto operator()(const nd::array_t<ProviderType, Rank>& value) const
    {
        return Dataspace::simple(shape(value));
    }
};

template<typename ProviderType, nd::uint Rank>
struct hdf5_datatype_creation<nd::array_t<ProviderType, Rank>>
{
    using value_type = std::decay_t<typename nd::array_t<ProviderType, Rank>::value_type>;

    auto operator()(const nd::array_t<ProviderType, Rank>& value) const
    {
        return make_datatype_for(value_type());
    }
};

template<typename ValueType, nd::uint Rank>
struct hdf5_container_creation<nd::shared_array<ValueType, Rank>>
{
    using value_type = ValueType;

    auto operator()(const Dataset& dset) const
    {
        auto fshape = dset.get_space().extent();
        auto mshape = nd::uivec_t<Rank>();

        if (size(fshape) != size(mshape))
            throw std::invalid_argument("HDF5 dataset has a different rank from the nd::array to be read");

        if (dset.get_type() != make_datatype_for(value_type()))
            throw std::invalid_argument("HDF5 dataset has a different type from the nd::array to be read");

        for (std::size_t i = 0; i < Rank; ++i)
            mshape[i] = fshape[i];

        return nd::make_unique_array<value_type, Rank>(mshape);
    }
};

template<typename ValueType, nd::uint Rank>
struct hdf5_container_conversion_post_read<nd::unique_array<ValueType, Rank>>
{
    auto operator()(nd::unique_array<ValueType, Rank>&& value) const
    {
        return make_array(std::move(value.provider).shared(), shape(value));
    }
};




//=============================================================================
template<typename T>
void write(const Group& group, std::string name, const std::array<nd::shared_array<T, 3>, 3>& v)
{
    write(group.require_group(name), "1", v.at(0));
    write(group.require_group(name), "2", v.at(1));
    write(group.require_group(name), "3", v.at(2));
}

} // namespace h5




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"
#include "core_numeric_tuple.hpp"
#include "core_numeric_array.hpp"
#include "core_dimensional.hpp"




//=============================================================================
inline void test_hdf5_ndarray()
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
    test_read_write(nd::linspace(0, 1, 100) | nd::to_shared());
    test_read_write(nd::range(1000) | nd::to_shared());
}

#endif // DO_UNIT_TESTS
