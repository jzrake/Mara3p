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
#include <cstring>




//=============================================================================
namespace serial {




/**
 * @brief      Unspecialized is-serializable class. In addition to specializing
 *             the appropriate type descriptor or conversion classes below, you
 *             must also specialize this class to inherit from std::true_type
 *             for the type you are serializing.
 *
 * @tparam     T     The type to serialize
 */
template<typename T>
struct is_serializable_t : std::false_type {};




/**
 * @brief      Unspecialized type-descriptor class. For each non-POD data type you
 *             want to serialize (types not in std::is_trivially_copyable) you
 *             must provide a template specializaton of this struct, which
 *             overloads the call operator as a template member function:
 *
 *
 *             template<typename Serializer>
 *             void operator()(Serializer& s, T& value) const
 *             {
 *                 s(value.member1);
 *                 s(value.member2);
 *             }
 *
 *
 * @tparam     T     The type to serialize
 */
template<typename T>
struct type_descriptor_t
{
};




/**
 * @brief      Unspecialized conversion-to-serializable class. You should
 *             specialize this struct for a type T if T must be converted to
 *             another type to be serialized. You must create a typedef called
 *             'type' to represent the converted-to type U in your
 *             specialization (obviously T must be a different type than U). You
 *             must also create an overload of the call operator taking a T and
 *             returning a U.
 *
 * @tparam     T     The serializable data type to convert to
 */
template<typename T>
struct conversion_to_serializable_t
{
    using type = T;
    auto operator()(T value) const { return value; }
};




/**
 * @brief      Unspecialized conversion-from-serializable class. You should
 *             specialize this struct for a type T if T is fetched from the
 *             deserializer as a type U other than T. You must create a typedef
 *             'type' to represent U, and overload the call operator to take a U
 *             and return a T.
 *
 * @tparam     T     The serializable data type that was converted to
 */
template<typename T>
struct conversion_from_serializable_t
{
    using type = T;
    auto operator()(T value) const { return value; }
};




/**
 * @brief      Unspecialized struct to write the size of a dynamically sized
 *             container type T. The call operator must be overloaded with the
 *             same signature as below, but where a description of the item size
 *             (or shape) is described to the Serializer instance using the same
 *             types of calls as in the type_descriptor_t's call operator. If
 *             this struct is specialized then you also need to specialize
 *             container_shape_setter_t to vend out the size or shape
 *             information that was described here.
 *
 * @tparam     T     The container type
 */
template<typename T>
struct container_shape_descriptor_t
{
    template<typename Serializer>
    void operator()(Serializer&, const T&) const {}
};




/**
 * @brief      Unspecialized struct to create an instance of a dynamically sized
 *             container type T. The call operator must be overloaded to vend
 *             from the Serializer instance the size or shape described in the
 *             specialization of container_shape_descriptor_t<T>, and then
 *             resize the input container appropriately.
 *
 * @tparam     T     The container type
 */
template<typename T>
struct container_shape_setter_t
{
    template<typename Serializer>
    void operator()(Serializer&, T&) const {}
};




//=============================================================================
template<typename T>
constexpr bool is_serializable()
{
    if constexpr (std::is_trivially_copyable_v<T>)
    {
        return true;
    }
    else
    {
        return is_serializable_t<T>::value;        
    }
}




//=============================================================================
class deserializer_t
{
public:

    deserializer_t(const std::vector<char>& buffer) : buffer(buffer) {}




    /**
     * @brief      Attempt to return an item of type T from the front of the
     *             internal buffer, and advance the position to point to the
     *             next serialized item.
     *
     * @tparam     T     The value type to vend out
     *
     * @return     An instance of T obtained from the internal buffer
     */
    template<typename T>
    T vend()
    {
        auto value = T();
        operator()(value);
        return value;
    }




    //=========================================================================
    template <typename T>
    void operator()(T& value)
    {
        if constexpr (std::is_trivially_copyable_v<T>)
        {
            std::memcpy(&value, buffer.data() + position, sizeof(T));
            position += sizeof(T);
        }
        else
        {
            do_deserialization(value);
        }
    }




    //=========================================================================
    template<typename T>
    void operator()(T* data, std::size_t count)
    {
        for (std::size_t i = 0; i < count; ++i)
        {
            operator()(data[i]);
        }
    }




private:
    //=========================================================================
    template<typename T>
    void do_deserialization(T& value)
    {
        if constexpr (! std::is_same_v<typename conversion_from_serializable_t<T>::type, T>)
        {
            value = conversion_from_serializable_t<T>()(vend<typename conversion_from_serializable_t<T>::type>());
        }
        else
        {
            container_shape_setter_t<T>()(*this, value);
            type_descriptor_t<T>()(*this, value);
        }
    }

    const std::vector<char>& buffer;
    std::size_t position = 0;
};




//=============================================================================
class serializer_t
{
public:




    //=========================================================================
    template<typename T>
    void operator()(T&& value)
    {
        operator()(value);
    }




    //=========================================================================
    template<typename T>
    void operator()(T& value)
    {
        if constexpr (std::is_trivially_copyable_v<T>)
        {
            buffer.resize(buffer.size() + sizeof(T));
            std::memcpy(buffer.data() + buffer.size() - sizeof(T), &value, sizeof(T));
        }
        else
        {
            container_shape_descriptor_t<T>()(*this, value);
            do_serialization(value);
        }
    }




    //=========================================================================
    template<typename T>
    void operator()(const T* data, std::size_t count)
    {
        for (std::size_t i = 0; i < count; ++i)
        {
            operator()(data[i]);
        }
    }




    /**
     * @brief      Get a const reference to the internal buffer
     *
     * @return     A reference to the buffer.
     */
    const std::vector<char>& get_buffer() const &
    {
        return buffer;
    }




    /**
     * @brief      Get the internal buffer, moved out of this object.
     *
     * @return     An rvalue reference to the buffer.
     */
    std::vector<char>&& get_buffer() &&
    {
        return std::move(buffer);
    }




private:
    //=========================================================================
    template<typename T>
    void do_serialization(T& value)
    {
        if constexpr (! std::is_same_v<typename conversion_to_serializable_t<T>::type, T>)
        {
            operator()(conversion_to_serializable_t<T>()(value));            
        }
        else
        {
            type_descriptor_t<T>()(*this, value);
        }
    }
    std::vector<char> buffer;
};




/**
 * @brief      Write a value to a std::vector<char> that can be de-serialized
 *             with the loads function below.
 *
 * @param[in]  value  The value to serialize
 *
 * @tparam     T      The value type
 *
 * @return     The buffer of serialized data
 */
template<typename T, typename = typename std::enable_if_t<is_serializable<T>()>>
auto dumps(T value)
{
    serializer_t ser;
    ser(value);
    return std::move(ser).get_buffer();
}




/**
 * @brief      Create a value by de-serializing it from the given buffer.
 *
 * @param[in]  buffer  The buffer, created with dumps
 *
 * @tparam     T       The type to de-serialize to
 *
 * @return     The de-serialized instance
 */
template<typename T, typename = typename std::enable_if_t<is_serializable<T>()>>
T loads(const std::vector<char>& buffer)
{
    return deserializer_t(buffer).vend<T>();
}

} // namespace serial




//=============================================================================
#ifdef DO_UNIT_TESTS
#include <vector>
#include <string>
#include "core_unit_test.hpp"




/**
 * @brief      Example support structs to read / write std::vector of any
 *             serializable type
 *
 * @tparam     T     The vector value type
 */
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
struct serial::is_serializable_t<std::vector<T>> : std::true_type {};

template<typename T>
struct serial::type_descriptor_t<std::vector<T>>
{
    template<typename Serializer>
    void operator()(Serializer& s, std::vector<T>& value) const
    {
        s(value.data(), value.size());
    }
};




/**
 * @brief      An example showing how to read / write strings by converting them
 *             to / from vector<char>
 */
template<>
struct serial::is_serializable_t<std::string> : std::true_type {};

template<>
struct serial::conversion_to_serializable_t<std::string>
{
    using type = std::vector<char>;
    auto operator()(std::string value) const { return std::vector<char>(value.begin(), value.end()); }
};

template<>
struct serial::conversion_from_serializable_t<std::string>
{
    using type = std::vector<char>;
    auto operator()(std::vector<char> value) const { return std::string(value.begin(), value.end()); }
};




/**
 * @brief      An example non-POD data structure that does not have
 *             serialization info
 */
struct not_serializable
{
    std::vector<int> values;
};





/**
 * @brief      An example POD data structure (does not need to specialize
 *             type_descriptor_t)
 */
struct pod_struct_t
{
    int a;
    double b;
};

inline bool operator==(pod_struct_t a, pod_struct_t b)
{
    return a.a == b.a && a.b == b.b;
}




/**
 * @brief      An example struct that is not POD, so it does specialize
 *             type_descriptor_t
 */
struct non_pod_struct_t
{
    int a;
    double b;
    std::vector<double> c;
    std::string d;
};

inline bool operator==(non_pod_struct_t a, non_pod_struct_t b)
{
    return a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d;
}

template<>
struct serial::is_serializable_t<non_pod_struct_t> : std::true_type {};

template<>
struct serial::type_descriptor_t<non_pod_struct_t>
{
    template<typename Serializer>
    void operator()(Serializer& s, non_pod_struct_t& value) const
    {
        s(value.a);
        s(value.b);
        s(value.c);
        s(value.d);
    }
};




//=============================================================================
inline void test_serial()
{
    auto require_serializes = [] (auto value)
    {
        require(value == serial::loads<decltype(value)>(serial::dumps(value)));
    };

    static_assert(! std::is_trivially_copyable_v<non_pod_struct_t>);
    require_serializes(std::vector<int>{0, 1, 2});
    require_serializes(std::string("hey there!"));
    require_serializes(pod_struct_t{32, 3.14159});
    require_serializes(non_pod_struct_t{32, 3.14159, {3.4, 0.0, 3.10101}, "does it work?"});

    require(  serial::is_serializable<int>());
    require(! serial::is_serializable<not_serializable>());
    require(  serial::is_serializable<non_pod_struct_t>());
    require(  serial::is_serializable<std::string>());
}

#endif // DO_UNIT_TESTS
