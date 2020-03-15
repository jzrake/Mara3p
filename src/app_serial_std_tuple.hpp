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
#include <array>
#include <tuple>
#include <utility>
#include "app_serial.hpp"




//=============================================================================
template<typename... T>
struct serial::type_descriptor_t<std::tuple<T...>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::tuple<T...>& value) const
    {
        std::apply([&s] (auto&&... vs) { (..., s(vs)); }, value);
    }
};

template<typename... T> struct serial::is_serializable_t<std::tuple<T...>> : std::true_type {};




//=============================================================================
template<typename T, typename U>
struct serial::type_descriptor_t<std::pair<T, U>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::pair<T, U>& value) const
    {
        std::apply([&s] (auto&&... vs) { (..., s(vs)); }, value);
    }
};

template<typename T, typename U> struct serial::is_serializable_t<std::pair<T, U>> : std::true_type {};




//=============================================================================
template<typename T, std::size_t S>
struct serial::type_descriptor_t<std::array<T, S>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::array<T, S>& value) const
    {
        std::apply([&s] (auto&&... vs) { (..., s(vs)); }, value);
    }
};

template<typename T, std::size_t S> struct serial::is_serializable_t<std::array<T, S>> : std::true_type {};
