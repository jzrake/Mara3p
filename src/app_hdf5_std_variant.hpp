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
#include <variant>
#include "app_hdf5.hpp"




//=============================================================================
namespace h5 {




//=============================================================================
template<typename... Ts>
struct hdf5_container_address<std::variant<Ts...>>
{
    const void* operator()(const std::variant<Ts...>& value) const
    {
        return std::visit([] (const auto& v) { return container_address_of(v); }, value);
    }
    void* operator()(std::variant<Ts...>& value) const
    {
        return std::visit([] (auto& v) { return container_address_of(v); }, value);
    }
};

template<typename... Ts>
struct hdf5_dataspace_creation<std::variant<Ts...>>
{
    auto operator()(const std::variant<Ts...>& value) const
    {
        return std::visit([] (auto v) { return make_dataspace_for(v); }, value);
    }
};

template<typename... Ts>
struct hdf5_datatype_creation<std::variant<Ts...>>
{
    auto operator()(const std::variant<Ts...>& value) const
    {
        return std::visit([] (auto v) { return make_datatype_for(v); }, value);
    }
};

template<typename... Ts>
struct hdf5_container_creation<std::variant<Ts...>>
{
    std::variant<Ts...> operator()(const Dataset& dset) const
    {
        auto matches = [&] (auto v) -> bool
        {
            if (dset.get_type().get_class() == H5T_STRING)
            {
                return dset.get_type().with_size(1) == make_datatype_for(v);
            }
            return dset.get_type() == make_datatype_for(v);
        };

        auto create = [&] (auto v) -> std::variant<Ts...>
        {
            return matches(v) ? create_container_for<decltype(v)>(dset) : v;
        };

        auto containers = std::apply([&] (auto... ts)
        {
            return std::array{std::pair(matches(ts), create(ts))...};
        }, std::tuple<Ts...>{});

        for (const auto& c : containers)
        {
            if (c.first)
            {
                return c.second;
            }
        }
        throw std::invalid_argument("h5::read<std::variant> (HDF5 type cannot be matched to the variant type)");
    }
};

} // namespace h5




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
void test_hdf5_std_variant()
{
    auto test_read_write = [] (auto value)
    {
        {
            auto file = h5::File("test.h5", h5::File::Access::Truncate);
            h5::write(file, "value", value);
        }
        {
            auto file = h5::File("test.h5", h5::File::Access::Read);
            require(h5::read<decltype(value)>(file, "value") == value);
        }
    };

    test_read_write(std::variant<int, double, std::string>{1.1});
    test_read_write(std::variant<int, double, std::string>{1});
    test_read_write(std::variant<int, double, std::string>{std::string("Hey there!")});
}

#endif // DO_UNIT_TESTS
