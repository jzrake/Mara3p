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
#include <vector>
#include "app_serial.hpp"




//=============================================================================
template<typename T>
struct serial::container_shape_setter_t<std::vector<T>>
{
    template<typename Serializer>
    auto operator()(Serializer& s, std::vector<T>& value)
    {
        value.resize(s.template vend<std::size_t>());
    }
};

template<typename T>
struct serial::container_shape_descriptor_t<std::vector<T>>
{
    template<typename Serializer>
    auto operator()(Serializer& s, const std::vector<T>& value)
    {
        s(value.size());
    }
};

template<typename T>
struct serial::type_descriptor_t<std::vector<T>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::vector<T>& value) const
    {
        s(value.data(), value.size());
    }
};


template<typename T>
struct serial::is_serializable_t<std::vector<T>> : std::true_type {};
