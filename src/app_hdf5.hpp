#include <optional>
#include <string>
#include <vector>
#include <hdf5.h>




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
        auto hdims = std::vector<hsize_t>(dims.begin(), dims.end());
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
        auto hdims = std::vector<hsize_t>(dims.begin(), dims.end());
        auto max_hdims = std::vector<hsize_t>(max_dims.begin(), max_dims.end());
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

    Datatype copy() const
    {
        return Identifier::datatype(H5Tcopy(id));
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
            H5Oget_info_by_name(id, name.data(), &info, H5P_DEFAULT);
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
            detail::check(H5Oget_info_by_name(id, name.data(), &info, H5P_DEFAULT));
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
     * @brief      Creates a dataset.
     *
     * @param[in]  name   The name
     * @param[in]  type   The type
     * @param[in]  space  The space
     *
     * @return     { description_of_the_return_value }
     */
    Dataset create_dataset(std::string name, const Datatype& type, const Dataspace& space) const
    {
        return Identifier::dataset(H5Dcreate(id, name.data(), type.id, space.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }




    /**
     * @brief      Opens a dataset.
     *
     * @param[in]  name  The name
     *
     * @return     { description_of_the_return_value }
     */
    Dataset open_dataset(std::string name) const
    {
        return Identifier::dataset(H5Dopen(id, name.data(), H5P_DEFAULT));
    }




    /**
     * @brief      { function_description }
     *
     * @param[in]  name   The name
     * @param[in]  type   The type
     * @param[in]  space  The space
     *
     * @return     { description_of_the_return_value }
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




    File() {}
    File(std::string filename, Access access)
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




/**
 * @brief      Unspecified conversion-to-serializable class. You should
 *             specialize this struct for a type T if T must be converted to
 *             another type to be serialized. You must create a typedef called
 *             'type' to represent the converted-to type U in your
 *             specialization (obviously T must be a different type than U). You
 *             must also create an overload of the call operator taking a T and
 *             returning a U.
 *
 * @tparam     T     The data type to convert from
 */
template<typename T>
struct conversion_to_serializable_t
{
    using type = T;
    auto operator()(T value) const { return value; }
};




/**
 * @brief      Unspecified conversion-from-serializable class. You should
 *             specialize this struct for a type T if T is fetched from the
 *             deserializer as a type U other than T. You must create a typedef
 *             'type' to represent U, and overload the call operator to take a U
 *             and return a T.
 *
 * @tparam     T     The data type to convert to
 */
template<typename T>
struct conversion_from_serializable_t
{
    using type = T;
    auto operator()(T value) const { return value; }
};




/**
 * @brief      Unspecified struct to write the size of a dynamically sized
 *             container type T. The call operator must be overloaded with the
 *             same signature as below, but where a description of the item size
 *             (or shape) is described to the Serializer instance using the same
 *             types of calls as in the type_descriptor_t's call operator. If
 *             this struct is specialized then you also need to specialize
 *             container_shape_setter_t to vend out the size or shape information
 *             that was described here.
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
 * @brief      Unspecified struct to create an instance of a dynamically sized
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
template<typename T> struct hdf5_container_address
{
    const void* operator()(const T& value) { return &value; }
    void* operator()(T& value) { return &value; }
};

template<typename T> struct hdf5_datatype_creation
{
};

template<typename T> struct hdf5_dataspace_creation
{
    Dataspace operator()(const T&) { return Dataspace::scalar(); }
};

template<typename T> struct hdf5_container_creation
{
    T operator()(const Dataset&) { return T(); }
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
T create_container_for(const Dataset& dataset)
{
    return hdf5_container_creation<T>()(dataset);
}




//=============================================================================
template<> struct hdf5_datatype_creation<double>        { Datatype operator()(const double&)        { return Datatype::native_double(); } };
template<> struct hdf5_datatype_creation<float>         { Datatype operator()(const float&)         { return Datatype::native_float(); } };
template<> struct hdf5_datatype_creation<int>           { Datatype operator()(const int&)           { return Datatype::native_int(); } };
template<> struct hdf5_datatype_creation<long>          { Datatype operator()(const long&)          { return Datatype::native_long(); } };
template<> struct hdf5_datatype_creation<unsigned long> { Datatype operator()(const unsigned long&) { return Datatype::native_ulong(); } };




//=============================================================================
template<> struct hdf5_datatype_creation<std::string>
{
    Datatype operator()(const std::string& v)
    {
        return Datatype::c_s1().with_size(std::max(std::size_t(1), v.size()));
    }
};

template<> struct hdf5_container_address<std::string>
{
    const void* operator()(const std::string& value) { return value.data(); }
    void* operator()(std::string& value) { return value.data(); }
};

template<> struct hdf5_container_creation<std::string>
{
    std::string operator()(const Dataset& dset) { return std::string(dset.get_type().size(), 0); }
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

template<typename T> struct hdf5_container_address<std::vector<T>>
{
    const void* operator()(const std::vector<T>& value) { return value.data(); }
    void* operator()(std::vector<T>& value) { return value.data(); }
};

template<typename T> struct hdf5_container_creation<std::vector<T>>
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
    auto type = make_datatype_for(value);
    auto mspc = make_dataspace_for(value);
    auto data = container_address_of(value);
    dset.write(type, mspc, mspc, data);
}

template<typename T>
void write(const Group& group, std::string name, const T& value)
{
    auto type = make_datatype_for(value);
    auto mspc = make_dataspace_for(value);
    write(group.require_dataset(name, type, mspc), value);
}

template<typename T>
void read(const Dataset& dset, T& value)
{
    value = create_container_for<T>(dset);
    auto type = make_datatype_for(value);
    auto fspc = dset.get_space();
    auto data = container_address_of(value);
    dset.read(type, fspc, fspc, data);
}

template<typename T>
void read(const Group& group, std::string name, T& value)
{
    read(group.open_dataset(name), value);
}

template<typename T>
T read(const Dataset& dset)
{
    auto value = T();
    read(dset, value);
    return value;
}

template<typename T>
T read(const Group& group, std::string name)
{
    return read<T>(group.open_dataset(name));
}

} // namespace h5




//=============================================================================
#ifdef DO_UNIT_TESTS
#include "core_unit_test.hpp"




//=============================================================================
void test_hdf5()
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
}

#endif // DO_UNIT_TESTS
