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
#include "app_serial.hpp"




//=============================================================================
namespace detail {

template<typename VariantType, std::size_t... I>
void construct_variant(std::size_t index, VariantType& value, std::index_sequence<I...>)
{
    (..., [&value, index] () { if (I == index) value = std::variant_alternative_t<I, VariantType>(); }());
}

}




//=============================================================================
template<typename... T>
struct serial::container_shape_descriptor_t<std::variant<T...>>
{
    template<typename Serializer>
    auto operator()(Serializer& s, const std::variant<T...>& value)
    {
        s(value.index());
    }
};

template<typename... T>
struct serial::container_shape_setter_t<std::variant<T...>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::variant<T...>& value) const
    {
        auto index = s.template vend<std::size_t>();
        detail::construct_variant<std::variant<T...>>(index, value, std::make_index_sequence<std::variant_size_v<std::variant<T...>>>());
    }
};

template<typename... T>
struct serial::type_descriptor_t<std::variant<T...>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::variant<T...>& value) const
    {
        std::visit([&s] (auto& v) { s(v); }, value);
    }
};

template<typename... T>
struct serial::is_serializable_t<std::variant<T...>> : std::true_type {};
