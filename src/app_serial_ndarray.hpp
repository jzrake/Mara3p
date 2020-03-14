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
#include "app_serial.hpp"
#include "core_ndarray.hpp"




//=============================================================================
template<typename ValueType, std::size_t Rank>
struct serial::conversion_from_serializable_t<nd::shared_array<ValueType, Rank>>
{
    using type = nd::unique_array<ValueType, Rank>;
    auto operator()(type value) const
    {
        return nd::make_array(std::move(value.provider).shared(), value.shape);
    }
};

template<typename ValueType, std::size_t Rank>
struct serial::container_shape_descriptor_t<nd::shared_array<ValueType, Rank>>
{
    template<typename Serializer>
    auto operator()(Serializer& s, const nd::shared_array<ValueType, Rank>& value)
    {
        s(shape(value));
    }
};

template<typename ValueType, std::size_t Rank>
struct serial::container_shape_setter_t<nd::unique_array<ValueType, Rank>>
{
    template<typename Serializer>
    void operator()(Serializer& s, nd::unique_array<ValueType, Rank>& value) const
    {
        value = nd::make_unique_array<ValueType>(s.template vend<nd::uivec_t<Rank>>());
    }
};

template<typename ValueType, std::size_t Rank>
struct serial::type_descriptor_t<nd::shared_array<ValueType, Rank>>
{
    template<typename Serializer>
    void operator()(Serializer& s, nd::shared_array<ValueType, Rank>& value) const
    {
        s(value.data(), size(value));
    }
};

template<typename ValueType, std::size_t Rank>
struct serial::type_descriptor_t<nd::unique_array<ValueType, Rank>>
{
    template<typename Serializer>
    void operator()(Serializer& s, nd::unique_array<ValueType, Rank>& value) const
    {
        s(value.data(), size(value));
    }
};

template<typename ValueType>
struct serial::type_descriptor_t<std::array<nd::shared_array<ValueType, 3>, 3>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::array<nd::shared_array<ValueType, 3>, 3>& value) const
    {
        s(value.at(0));
        s(value.at(1));
        s(value.at(2));
    }
};




//=============================================================================
template<typename ValueType, std::size_t Rank> struct serial::is_serializable_t<nd::shared_array<ValueType, Rank>> : std::true_type {};
template<typename ValueType, std::size_t Rank> struct serial::is_serializable_t<nd::unique_array<ValueType, Rank>> : std::true_type {};
template<typename ValueType> struct serial::is_serializable_t<std::array<nd::shared_array<ValueType, 3>, 3>> : std::true_type {};




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "app_serial_numeric_tuple.hpp"
#include "core_dimensional.hpp"
#include "core_numeric_array.hpp"
#include "core_numeric_tuple.hpp"
#include "core_unit_test.hpp"




//=============================================================================
inline void test_serial_ndarray()
{
    auto test_serializes = [] (auto v1)
    {
        auto buffer = serial::dumps(v1);
        require(all(serial::loads<decltype(v1)>(buffer) == v1));
    };
    auto test_serializes2 = [] (auto v1)
    {
        auto buffer = serial::dumps(v1);
        require(serial::loads<decltype(v1)>(buffer) == v1);
    };
    test_serializes(nd::linspace(0, 1, 100) | nd::to_shared());
    test_serializes(nd::range(1000) | nd::to_shared());
    test_serializes(nd::from(numeric::tuple(1, 2, 3)) | nd::to_shared());
    test_serializes(nd::from(numeric::array(1, 2, 3)) | nd::to_shared());
    test_serializes(nd::from(numeric::array(dimensional::unit_mass(1), dimensional::unit_mass(2), dimensional::unit_mass(3))) | nd::to_shared());
    test_serializes(nd::from(numeric::tuple(2, 1.0)) | nd::to_shared());
    test_serializes2(numeric::tuple(1, 2.0, 3));
}

#endif // DO_UNIT_TESTS
