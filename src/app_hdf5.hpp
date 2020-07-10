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
#include <tuple>
#include <optional>
#include <string>
#include <vector>
#include <hdf5.h>




//=============================================================================
#define h5_compound_type_member(type, name) {#name, offsetof(type, name), h5::make_datatype_for(type().name)}




//=============================================================================
namespace h5 {

namespace detail {

inline herr_t get_last_error(unsigned n, const H5E_error2_t *err, void *data)
{
    if (n == 0)
    {
        *static_cast<H5E_error2_t*>(data) = *err;
    }
    return 0;
}

inline hid_t check(hid_t result)
{
    if (result < 0)
    {
        H5E_error2_t err;
        hid_t eid = H5Eget_current_stack();
        H5Ewalk(eid, H5E_WALK_UPWARD, get_last_error, &err);
        std::string what = err.desc;
        H5Eclear(eid);
        H5Eclose_stack(eid);
        throw std::invalid_argument(what);
    }
    return result;
}

}




//=============================================================================
class Identifier
{
public:




    /**
     * @brief      This class describes one of the HDF5 identifier types.
     */
    enum class Type
    {
        Dataset,
        Dataspace,
        Datatype,
        File,
        Group,
    };




    Identifier() {}
    Identifier(const Identifier& other) = delete;
    Identifier(Identifier&& other) : id(other.id)
    {
        other.id.reset();
    }
    ~Identifier() { close(); }




    /**
     * @brief      Move-assignment operator for an identifier.
     *
     * @param      other  The identifier to be moved from
     *
     * @return     This identifier after the assignment
     */
    Identifier& operator=(Identifier&& other)
    {
        close();
        id = other.id;
        other.id.reset();
        return *this; 
    }




    /**
     * @brief      Close this identifier, decrementing its reference count
     *             within the HDF5 library.
     */
    void close()
    {
        if (id.has_value())
        {
            switch (id.value().second)
            {
                case Type::Dataset:    H5Dclose(id.value().first); break;
                case Type::Dataspace:  H5Sclose(id.value().first); break;
                case Type::Datatype:   H5Tclose(id.value().first); break;
                case Type::File:       H5Fclose(id.value().first); break;
                case Type::Group:      H5Gclose(id.value().first); break;
            }
            id.reset();
        }
    }




    /**
     * @brief      Determine whether this identifier is open.
     *
     * @return     True if open, False otherwise.
     */
    bool is_open() const
    {
        return id.has_value();
    }



private:



    /**
     * @brief      Get the untyped HDF5 identifier, if this is open. Otherwise
     *             throws an error.
     *
     * @return     The HDF5 id.
     */
    operator hid_t() const
    {
        if (! is_open())
        {
            throw std::invalid_argument("h5::Identifier (id is not open)");
        }
        return id.value().first;
    }




    //=========================================================================
    friend class Dataset;
    friend class Dataspace;
    friend class Datatype;
    friend class File;
    friend class Group;

    Identifier(hid_t id, Type type) : id(std::pair(detail::check(id), type)) {}
    static auto dataset  (hid_t id) { return Identifier(id, Type::Dataset); }
    static auto dataspace(hid_t id) { return Identifier(id, Type::Dataspace); }
    static auto datatype (hid_t id) { return Identifier(id, Type::Datatype); }
    static auto file     (hid_t id) { return Identifier(id, Type::File); }
    static auto group    (hid_t id) { return Identifier(id, Type::Group); }

    std::optional<std::pair<hid_t, Type>> id;
};




//=============================================================================
class Dataspace
{
public:




    /**
     * @brief      Return a new scalar dataspace
     *
     * @return     A new dataspace
     */
    static Dataspace scalar()
    {
        return Identifier::dataspace(H5Screate(H5S_SCALAR));
    }




    /**
     * @brief      Return a new, non-resizable, n-dimensional dataspace from a
     *             container.
     *
     * @param[in]  dims       The dimensions
     *
     * @tparam     Container  The container type
     *
     * @return     A new dataspace
     */
    template<typename Container>
    static Dataspace simple(Container dims)
    {
        auto hdims = std::vector<hsize_t>(begin(dims), end(dims));
        return Identifier::dataspace(H5Screate_simple(hdims.size(), &hdims[0], nullptr));
    }




    /**
     * @brief      Return a new, resizable, scalar dataspace with the given
     *             maximum size
     *
     * @param[in]  dims       The dimensions
     * @param[in]  max_dims   The maximum dimensions
     *
     * @tparam     Container  The container type
     *
     * @return     A new dataspace
     */
    template<typename Container>
    static Dataspace simple(Container dims, Container max_dims)
    {
        if (dims.size() != max_dims.size())
        {
            throw std::invalid_argument("h5::Dataspace::simple (dims and max dims sizes do not agree)");
        }
        auto hdims = std::vector<hsize_t>(begin(dims), end(dims));
        auto max_hdims = std::vector<hsize_t>(begin(max_dims), end(max_dims));
        return Identifier::dataspace(H5Screate_simple(hdims.size(), &hdims[0], &max_hdims[0]));
    }




    /**
     * @brief      Return a new dataspace with unlimited size on all axes, and
     *             the given initial sizes.
     *
     * @param[in]  initial_sizes  The initial sizes for the dataspace
     *
     * @tparam     Args           The argument types
     *
     * @return     A new dataspace
     */
    template<typename... Args>
    static Dataspace unlimited(Args... initial_sizes)
    {
        auto initial = std::vector<hsize_t> {hsize_t(initial_sizes)...};
        auto maximum = std::vector<hsize_t>(sizeof...(Args), H5S_UNLIMITED);
        return simple(initial, maximum);
    }




    Dataspace() {}
    Dataspace(Identifier&& id) : id(std::move(id)) {}


    bool operator==(const Dataspace& other) const { return   detail::check(H5Sextent_equal(id, other.id)); }
    bool operator!=(const Dataspace& other) const { return ! detail::check(H5Sextent_equal(id, other.id)); }




    /**
     * @brief      Return the rank of this dataspace
     *
     * @return     The rank
     */
    std::size_t rank() const
    {
        return detail::check(H5Sget_simple_extent_ndims(id));
    }




    /**
     * @brief      Return the number of points in this dataspace
     *
     * @return     The number of points
     */
    std::size_t size() const
    {
        return detail::check(H5Sget_simple_extent_npoints(id));
    }




    /**
     * @brief      Return the extent (shape) of this dataspace
     *
     * @return     The dataspace extent
     */
    std::vector<std::size_t> extent() const
    {
        auto ext = std::vector<hsize_t>(rank());
        detail::check(H5Sget_simple_extent_dims(id, &ext[0], nullptr));
        return std::vector<std::size_t>(ext.begin(), ext.end());
    }




private:
    friend class Dataset;
    friend class Group;
    Identifier id;
};




//=============================================================================
class Datatype
{
public:


    //=========================================================================
    static Datatype native_double() { return Identifier(H5Tcopy(H5T_NATIVE_DOUBLE), Identifier::Type::Datatype); }
    static Datatype native_float()  { return Identifier(H5Tcopy(H5T_NATIVE_FLOAT),  Identifier::Type::Datatype); }
    static Datatype native_int()    { return Identifier(H5Tcopy(H5T_NATIVE_INT),    Identifier::Type::Datatype); }
    static Datatype native_long()   { return Identifier(H5Tcopy(H5T_NATIVE_LONG),   Identifier::Type::Datatype); }
    static Datatype native_ulong()  { return Identifier(H5Tcopy(H5T_NATIVE_ULONG),  Identifier::Type::Datatype); }
    static Datatype c_s1()          { return Identifier(H5Tcopy(H5T_C_S1),          Identifier::Type::Datatype); }


    //=========================================================================
    template<typename StructureType>
    static Datatype compound(std::initializer_list<std::tuple<const char*, std::ptrdiff_t, Datatype>> members)
    {
        static_assert(std::is_standard_layout<StructureType>::value,
            "can only make compound HDF5 types from standard layout structs");

        auto h5_type = detail::check(H5Tcreate(H5T_COMPOUND, sizeof(StructureType)));

        for (const auto& member : members)
        {
            detail::check(H5Tinsert(h5_type, std::get<0>(member), std::get<1>(member), std::get<2>(member).id));
        }
        return Identifier::datatype(h5_type);
    }


    Datatype() {}
    Datatype(Identifier&& id) : id(std::move(id)) {}


    bool operator==(const Datatype& other) const { return   detail::check(H5Tequal(id, other.id)); }
    bool operator!=(const Datatype& other) const { return ! detail::check(H5Tequal(id, other.id)); }


    std::size_t size() const
    {
        return detail::check(H5Tget_size(id));
    }

    bool is_variable_string() const
    {
        return detail::check(H5Tis_variable_str(id));
    }

    auto get_class() const
    {
        return detail::check(H5Tget_class(id));
    }

    Datatype copy() const
    {
        return Identifier::datatype(H5Tcopy(id));
    }

    Datatype get_native_type() const
    {
        return Identifier::datatype(H5Tget_native_type(id, H5T_DIR_ASCEND));
    }

    Datatype with_variable_size() const
    {
        auto result = copy();
        detail::check(H5Tset_size(result.id, H5T_VARIABLE));
        return result;        
    }

    Datatype with_size(std::size_t size) const
    {
        auto result = copy();
        detail::check(H5Tset_size(result.id, size));
        return result;
    }

    Datatype as_array(std::size_t size) const
    {
        auto dims = hsize_t(size);
        return Identifier::datatype(H5Tarray_create(id, 1, &dims));
    }

private:
    friend class Dataset;
    friend class Group;
    Identifier id;
};




//=============================================================================
class Dataset
{
public:
    Dataset() {}
    Dataset(Identifier&& id) : id(std::move(id)) {}

    Dataspace get_space() const
    {
        return Identifier::dataspace(H5Dget_space(id));
    }

    Datatype get_type() const
    {
        return Identifier::datatype(H5Dget_type(id));
    }

    void write(const Datatype& type, const Dataspace& mspace, const Dataspace& fspace, const void* data) const
    {
        detail::check(H5Dwrite(id, type.id, mspace.id, fspace.id, H5P_DEFAULT, data));
    }

    void read(const Datatype& type, const Dataspace& mspace, const Dataspace& fspace, void* data) const
    {
        detail::check(H5Dread(id, type.id, mspace.id, fspace.id, H5P_DEFAULT, data));
    }

private:
    Identifier id;
};




//=============================================================================
class Group
{
public:




    Group() {}
    Group(Identifier&& id) : id(std::move(id)) {}




    /**
     * @brief      Return true if this group has a subgroup with the given name
     *
     * @param[in]  name  The name of the group
     *
     * @return     True or false
     */
    bool contains_group(std::string name) const
    {
        if (H5Lexists(id, name.data(), H5P_DEFAULT))
        {
            H5O_info_t info;
#if H5_VERS_MINOR >= 12
            detail::check(H5Oget_info_by_name(id, name.data(), &info, 1, H5P_DEFAULT));
#else
            detail::check(H5Oget_info_by_name(id, name.data(), &info, H5P_DEFAULT));
#endif
            return info.type == H5O_TYPE_GROUP;
        }
        return false;
    }




    /**
     * @brief      Return true if this group has a dataset with the given name
     *
     * @param[in]  name  The name of the dataset
     *
     * @return     True or false
     */
    bool contains_dataset(std::string name) const
    {
        if (H5Lexists(id, name.data(), H5P_DEFAULT))
        {
            H5O_info_t info;
#if H5_VERS_MINOR >= 12
            detail::check(H5Oget_info_by_name(id, name.data(), &info, 1, H5P_DEFAULT));
#else
            detail::check(H5Oget_info_by_name(id, name.data(), &info, H5P_DEFAULT));
#endif
            return info.type == H5O_TYPE_DATASET;
        }
        return false;
    }




    /**
     * @brief      Create a group with the given name. Throws an exception if
     *             the group already exists.
     *
     * @param[in]  group_name  The group name to open
     *
     * @return     The group
     */
    Group create_group(std::string group_name) const
    {
        return Identifier::group(H5Gcreate(id, group_name.data(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }




    /**
     * @brief      Open a group with the given name. Throws an exception if the
     *             group does not exist.
     *
     * @param[in]  group_name  The group name to open
     *
     * @return     The group
     */
    Group open_group(std::string group_name) const
    {
        return Identifier::group(H5Gopen(id, group_name.data(), H5P_DEFAULT));
    }




    /**
     * @brief      Open or create a group with the given name.
     *
     * @param[in]  group_name  The group name to open or create
     *
     * @return     The group
     */
    Group require_group(std::string group_name) const
    {
        return contains_group(group_name) ? open_group(group_name) : create_group(group_name);
    }




    /**
     * @brief      Creates a dataset. Throws an exception if a dataset with the
     *             given name already exists.
     *
     * @param[in]  name   The name of the new dataset
     * @param[in]  type   The type of the new dataset
     * @param[in]  space  The shape of the new dataset
     *
     * @return     A dataset
     */
    Dataset create_dataset(std::string name, const Datatype& type, const Dataspace& space) const
    {
        return Identifier::dataset(H5Dcreate(id, name.data(), type.id, space.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }




    /**
     * @brief      Opens a dataset that already exists.
     *
     * @param[in]  name  The name of the dataset
     *
     * @return     A dataset
     */
    Dataset open_dataset(std::string name) const
    {
        return Identifier::dataset(H5Dopen(id, name.data(), H5P_DEFAULT));
    }




    /**
     * @brief      Either open or create a dataset with the given name, type,
     *             and shape. Throws an exception if an incompatible dataset
     *             with that name already exists.
     *
     * @param[in]  name   The name of the dataset to open or create
     * @param[in]  type   The type of the dataset
     * @param[in]  space  The shape of the dataset
     *
     * @return     A dataset
     */
    Dataset require_dataset(std::string name, const Datatype& type, const Dataspace& space) const
    {
        if (contains_dataset(name))
        {
            if (auto dset = open_dataset(name); dset.get_type() == type && dset.get_space() == space)
            {
                return dset;
            }
            throw std::invalid_argument("h5::Group::require_dataset (incompatible dataset already exists)");
        }
        return create_dataset(name, type, space);
    }




    //=========================================================================
    class iterator
    {
    public:

        using value_type = std::string;
        using iterator_category = std::forward_iterator_tag;

        //=====================================================================
        iterator(hid_t id, hsize_t idx) : id(id), idx(idx) {}
        iterator& operator++() { ++idx; return *this; }
        iterator operator++(int) { auto ret = *this; this->operator++(); return ret; }
        bool operator==(iterator other) const { return id == other.id && idx == other.idx; }
        bool operator!=(iterator other) const { return id != other.id || idx != other.idx; }
        std::string operator*() const
        {
            static char name[1024];

            if (H5Lget_name_by_idx(id, ".",
                H5_INDEX_NAME, H5_ITER_NATIVE, idx, name, 1024, H5P_DEFAULT) > 1024)
            {
                throw std::overflow_error("object names longer than 1024 are not supported");
            }
            return name;
        }

    private:
        //=====================================================================
        hid_t id = -1;
        hsize_t idx = 0;
    };




    /**
     * @brief      Return the number of items in this group.
     *
     * @return     The number
     */
    std::size_t size() const
    {
        auto op = [] (auto, auto, auto, auto) { return 0; };
        auto idx = hsize_t(0);
        H5Literate(id, H5_INDEX_NAME, H5_ITER_INC, &idx, op, nullptr);
        return idx;
    }

    iterator begin() const { return iterator(id, 0); }
    iterator end()   const { return iterator(id, size()); }




private:
    Identifier id;
};




//=============================================================================
class File
{
public:




    /**
     * @brief      This class describes an HDF5 access mode.
     */
    enum class Access
    {
        Read,
        ReadWrite,
        Truncate,
    };



    static Access access(std::string mode)
    {
        if (mode == "r")  return Access::Read;
        if (mode == "r+") return Access::ReadWrite;
        if (mode == "w")  return Access::Truncate;

        throw std::invalid_argument("h5::File (no such access mode " + mode + ")");
    }



    File() {}
    File(std::string filename, std::string mode) : File(filename, access(mode)) {}
    File(std::string filename, Access access=Access::Read)
    {
        switch (access)
        {
            case Access::Read:      id = Identifier::file(H5Fopen  (filename.data(), H5F_ACC_RDONLY, H5P_DEFAULT)); break;
            case Access::ReadWrite: id = Identifier::file(H5Fopen  (filename.data(), H5F_ACC_RDWR,   H5P_DEFAULT)); break;
            case Access::Truncate:  id = Identifier::file(H5Fcreate(filename.data(), H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT)); break;
        }
    }




    /**
     * @brief      Close the file handle.
     */
    void close()
    {
        id = {};
    }




    /**
     * @brief      Convert this file to its root group.
     */
    operator Group() const
    {
        return Identifier::group(H5Gopen(id, "/", H5P_DEFAULT));
    }




private:
    Identifier id;
};

} // namespace h5




//=============================================================================
namespace h5 {




//=============================================================================
template<typename T>
struct hdf5_container_address
{
    const void* operator()(const T& value) const { return &value; }
    void* operator()(T& value) const { return &value; }
};

template<typename T>
struct hdf5_datatype_creation
{
};

template<typename T>
struct hdf5_dataspace_creation
{
    Dataspace operator()(const T&) const { return Dataspace::scalar(); }
};

template<typename T>
struct hdf5_container_creation
{
    T operator()(const Dataset&) const { return T(); }
};

template<typename T>
struct hdf5_container_conversion_post_read
{
    auto operator()(T&& value) const { return std::move(value); }
};

template<typename T>
struct hdf5_is_key_value_container : std::false_type
{
};

template<typename T>
struct hdf5_conversion_to_hdf5_writable
{
    using type = T;
};

template<typename T>
struct hdf5_conversion_from_hdf5_writable
{
    using type = T;
};




//=============================================================================
template<typename T>
Datatype make_datatype_for(const T& value)
{
    return hdf5_datatype_creation<T>()(value);
}

template<typename T>
Dataspace make_dataspace_for(const T& value)
{
    return hdf5_dataspace_creation<T>()(value);
}

template<typename T>
const void* container_address_of(const T& value)
{
    return hdf5_container_address<T>()(value);
}

template<typename T>
void* container_address_of(T& value)
{
    return hdf5_container_address<T>()(value);
}

template<typename T>
auto create_container_for(const Dataset& dataset)
{
    return hdf5_container_creation<T>()(dataset);
}

template<typename T>
auto convert_container_post_read(T&& value)
{
    return hdf5_container_conversion_post_read<T>()(std::move(value));
}

template<typename T>
auto convert_to_hdf5_writable(const T& value)
{
    return hdf5_conversion_to_hdf5_writable<T>()(value);
}

template<typename T>
auto convert_from_hdf5_writable(const typename hdf5_conversion_from_hdf5_writable<T>::type& value)
{
    return hdf5_conversion_from_hdf5_writable<T>()(value);
}




//=============================================================================
template<> struct hdf5_datatype_creation<double>        { Datatype operator()(const double&)        { return Datatype::native_double(); } };
template<> struct hdf5_datatype_creation<float>         { Datatype operator()(const float&)         { return Datatype::native_float(); } };
template<> struct hdf5_datatype_creation<int>           { Datatype operator()(const int&)           { return Datatype::native_int(); } };
template<> struct hdf5_datatype_creation<long>          { Datatype operator()(const long&)          { return Datatype::native_long(); } };
template<> struct hdf5_datatype_creation<unsigned long> { Datatype operator()(const unsigned long&) { return Datatype::native_ulong(); } };




//=============================================================================
template<>
struct hdf5_datatype_creation<std::string>
{
    Datatype operator()(const std::string& v)
    {
        return Datatype::c_s1().with_size(std::max(std::size_t(1), v.size()));
    }
};

template<>
struct hdf5_container_address<std::string>
{
    const void* operator()(const std::string& value) { return value.data(); }
    void* operator()(std::string& value) { return value.data(); }
};

template<>
struct hdf5_container_creation<std::string>
{
    std::string operator()(const Dataset& dset)
    {
        return std::string(dset.get_type().size(), 0);
    }
};




//=============================================================================
template<typename T>
struct hdf5_dataspace_creation<std::vector<T>>
{
    auto operator()(const std::vector<T>& value) const
    {
        return Dataspace::simple(std::vector{value.size()});
    }
};

template<typename T>
struct hdf5_datatype_creation<std::vector<T>>
{
    auto operator()(const std::vector<T>& value) const
    {
        return make_datatype_for<T>({});
    }
};

template<typename T>
struct hdf5_container_address<std::vector<T>>
{
    const void* operator()(const std::vector<T>& value) { return value.data(); }
    void* operator()(std::vector<T>& value) { return value.data(); }
};

template<typename T>
struct hdf5_container_creation<std::vector<T>>
{
    std::vector<T> operator()(const Dataset& dset)
    {
        return std::vector<T>(dset.get_space().size());
    }
};




//=============================================================================
template<typename T>
void write(const Dataset& dset, const T& value)
{
    if constexpr (! std::is_same_v<typename hdf5_conversion_to_hdf5_writable<T>::type, T>)
    {
        write(dset, convert_to_hdf5_writable(value));
    }
    else
    {
        static_assert(! hdf5_is_key_value_container<T>::value);

        auto type = make_datatype_for(value);
        auto mspc = make_dataspace_for(value);
        auto data = container_address_of(value);
        dset.write(type, mspc, mspc, data);
    }
}

template<typename T>
void write(const Group& group, std::string name, const T& value)
{
    if constexpr (! std::is_same_v<typename hdf5_conversion_to_hdf5_writable<T>::type, T>)
    {
        write(group, name, convert_to_hdf5_writable(value));
    }
    else if constexpr (hdf5_is_key_value_container<T>::value)
    {
        auto subgroup = group.require_group(name);

        for (const auto& item : value)
        {
            write(subgroup, item.first, item.second);
        }
    }
    else
    {
        auto type = make_datatype_for(value);
        auto mspc = make_dataspace_for(value);
        write(group.require_dataset(name, type, mspc), value);
    }
}

template<typename T>
void read(const Dataset& dset, T& value)
{
    if constexpr (! std::is_same_v<typename hdf5_conversion_from_hdf5_writable<T>::type, T>)
    {
        auto value_in_file = typename hdf5_conversion_from_hdf5_writable<T>::type();
        read(dset, value_in_file);
        value = convert_from_hdf5_writable<T>(value_in_file);
    }
    else
    {
        static_assert(! hdf5_is_key_value_container<T>::value);

        auto temp = create_container_for<T>(dset);
        auto type = make_datatype_for(temp);
        auto fspc = dset.get_space();
        auto data = container_address_of(temp);
        dset.read(type, fspc, fspc, data);
        value = convert_container_post_read(std::move(temp));
    }
}

template<typename T>
void read(const Group& group, std::string name, T& value)
{
    if constexpr (hdf5_is_key_value_container<T>::value)
    {
        auto subgroup = group.open_group(name);

        for (auto key : subgroup)
        {
            value[key] = read<typename T::mapped_type>(subgroup, key);
        }
    }
    else
    {
        read(group.open_dataset(name), value);
    }
}

template<typename T>
T read(const Dataset& dset)
{
    static_assert(! hdf5_is_key_value_container<T>::value);

    auto value = T();
    read(dset, value);
    return value;
}

template<typename T>
T read(const Group& group, std::string name)
{
    if constexpr (hdf5_is_key_value_container<T>::value)
    {
        auto value = T();
        read(group, name, value);
        return value;
    }
    else
    {
        return read<T>(group.open_dataset(name));
    }
}

} // namespace h5




//=============================================================================
#ifdef DO_UNIT_TESTS
#include <tuple>
#include "core_unit_test.hpp"




template<>
struct h5::hdf5_conversion_to_hdf5_writable<std::tuple<int>>
{
    using type = int;
    auto operator()(const std::tuple<int>& value) const { return std::get<0>(value); }
};

template<>
struct h5::hdf5_conversion_from_hdf5_writable<std::tuple<int>>
{
    using type = int;
    auto operator()(const int& value) const { return std::tuple(value); }
};




//=============================================================================
inline void test_hdf5()
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

    test_read_write(10.0);
    test_read_write(10.0f);
    test_read_write(10);
    test_read_write(10l);
    test_read_write(10ul);
    test_read_write(std::string("ten"));
    test_read_write(std::vector{1, 2, 3});
    test_read_write(std::vector{1.1, 2.2, 3.3});
    test_read_write(std::tuple(13));

    {
        auto file = h5::File("test.h5", h5::File::Access::Truncate);
        h5::Group(file).create_group("group");
        require(  h5::Group(file).contains_group("group"));
        require(! h5::Group(file).contains_group("froup"));
    }
}

#endif // DO_UNIT_TESTS
