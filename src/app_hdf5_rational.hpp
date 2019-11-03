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
#include "core_rational.hpp"




//=============================================================================
template<>
struct h5::hdf5_datatype_creation<rational::number_t>
{
    auto operator()(const rational::number_t& value) const
    {
        return Datatype::compound<rational::number_t>({
            h5_compound_type_member(rational::number_t, num),
            h5_compound_type_member(rational::number_t, den),
        });
    }
};




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "app_hdf5_rational.hpp"
#include "core_unit_test.hpp"




//=============================================================================
inline void test_hdf5_rational()
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
    test_read_write(rational::number(1));
    test_read_write(rational::number(3, 2));
}

#endif // DO_UNIT_TESTS
