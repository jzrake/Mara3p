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




#include <ostream>
#include "core_dimensional.hpp"
#include "core_geometric.hpp"
#include "core_ndarray.hpp"




//=============================================================================
namespace vtk {




//=============================================================================
using quadratic_cell_t = std::array<int, 4>;




//=============================================================================
template<typename T, typename = std::enable_if_t<std::is_trivially_copyable_v<T>>>
std::array<char, sizeof(T)> swap_bytes(T value)
{
    auto result = std::array<char, sizeof(T)>();

    for (std::size_t i = 0; i < sizeof(T); ++i)
        result[sizeof(T) - i - 1] = reinterpret_cast<char*>(&value)[i];

    return result;
}

template<typename T, typename = std::enable_if_t<std::is_trivially_copyable_v<T>>>
void write_swapped(std::ostream& os, T value)
{
    os.write(swap_bytes(value).data(), sizeof(T));
}

template<typename P, typename FunctionType>
void column_major_loop(nd::array_t<P, 3> A, FunctionType f)
{
    for (nd::uint k = 0; k < shape(A, 2); ++k)
        for (nd::uint j = 0; j < shape(A, 1); ++j)
            for (nd::uint i = 0; i < shape(A, 0); ++i)
                f(A(i, j, k));
}




//=============================================================================
void write_structured_grid_header(std::ostream& os, std::string title, nd::uivec_t<3> vertices_shape)
{
    auto [ni, nj, nk] = vertices_shape.impl;
    os << "# vtk DataFile Version 3.0\n";
    os << title << '\n';
    os << "BINARY\n";
    os << "DATASET STRUCTURED_GRID\n";
    os << "DIMENSIONS " << ni << ' ' << nj << ' ' << nk << '\n';
}

void write_unstructured_grid_header(std::ostream& os, std::string title, nd::uint cell_count)
{
    os << "# vtk DataFile Version 3.0\n";
    os << title << '\n';
    os << "BINARY\n";
    os << "DATASET UNSTRUCTURED_GRID\n";
}

void write_point_data_prelude(std::ostream& os, nd::uivec_t<3> points_shape)
{
    os << "POINT_DATA " << product(points_shape) << '\n';
}

void write_cell_data_prelude(std::ostream& os, nd::uivec_t<3> cells_shape)
{
    os << "CELL_DATA " << product(cells_shape) << '\n';
}

void write_cell_data_prelude(std::ostream& os, nd::uint cell_count)
{
    os << "CELL_DATA " << cell_count << '\n';
}




//=============================================================================
template<typename T, typename = std::enable_if_t<sizeof(T) == sizeof(double)>>
void write_point_coordinates(std::ostream& os, nd::shared_array<geometric::euclidean_vector_t<T>, 3> X)
{
    os << "POINTS " << size(X) << " double\n";

    column_major_loop(X, [&os] (const auto& x)
    {
        write_swapped(os, x.component_1());
        write_swapped(os, x.component_2());
        write_swapped(os, x.component_3());
    });
}

template<typename T, typename = std::enable_if_t<sizeof(T) == sizeof(double)>>
inline void write_point_coordinates(std::ostream& os, nd::shared_array<geometric::euclidean_vector_t<T>, 1> X)
{
    os << "POINTS " << size(X) << " double\n";

    for (auto x : X)
    {
        write_swapped(os, x.component_1());
        write_swapped(os, x.component_2());
        write_swapped(os, x.component_3());
    }
}

inline void write_cells(std::ostream& os, nd::shared_array<quadratic_cell_t, 1> X)
{
    int num_points_per_cell = 4;
    int quadratic_cell_type = 9;

    os << "CELLS " << size(X) << " " << size(X) * 5 << '\n';

    for (auto q : X)
    {
        write_swapped(os, num_points_per_cell);
        write_swapped(os, q[0]);
        write_swapped(os, q[1]);
        write_swapped(os, q[2]);
        write_swapped(os, q[3]);
    }

    os << "CELL_TYPES " << size(X) << '\n';

    for (auto q : X)
    {
        (void)q;
        write_swapped(os, quadratic_cell_type);
    }
}




//=============================================================================
template<typename T, typename = std::enable_if_t<sizeof(T) == sizeof(double)>>
void write_cell_data(std::ostream& os, std::string name, nd::shared_array<T, 3> A)
{
    os << "SCALARS " << name << " double\n";
    os << "LOOKUP_TABLE default\n";

    column_major_loop(A, [&os] (const auto& a)
    {
        write_swapped(os, a);
    });
}

//=============================================================================
template<typename T, typename = std::enable_if_t<sizeof(T) == sizeof(double)>>
void write_cell_data(std::ostream& os, std::string name, nd::shared_array<T, 1> A)
{
    os << "SCALARS " << name << " double\n";
    os << "LOOKUP_TABLE default\n";

    for (auto a : A)
    {
        write_swapped(os, a);
    }
}

//=============================================================================
template<typename T, typename = std::enable_if_t<sizeof(T) == sizeof(double)>>
void write_cell_data(std::ostream& os, std::string name, nd::shared_array<geometric::euclidean_vector_t<T>, 3> V)
{
    os << "VECTORS " << name << " double\n";

    column_major_loop(V, [&os] (const auto& v)
    {
        write_swapped(os, v.component_1());
        write_swapped(os, v.component_2());
        write_swapped(os, v.component_3());
    });
}




//=============================================================================
template<typename VertexPositionType, typename... CellDataArrayTypes>
void write(std::ostream& os,
    const char* dataset_name,
    nd::shared_array<geometric::euclidean_vector_t<VertexPositionType>, 3> vertices,
    std::pair<const char*, CellDataArrayTypes>... cell_fields)
{
    auto [ni, nj, nk] = vertices.shape.impl;

    write_structured_grid_header(os, dataset_name, shape(vertices));
    write_point_coordinates(os, vertices);
    write_cell_data_prelude(os, {ni - 1, nj - 1, nk - 1});
    (..., write_cell_data(os, cell_fields.first, cell_fields.second));
}




//=============================================================================
template<typename VertexPositionType, typename... CellDataArrayTypes>
void write(std::ostream& os,
    const char* dataset_name,
    nd::shared_array<geometric::euclidean_vector_t<VertexPositionType>, 1> vertices,
    nd::shared_array<quadratic_cell_t, 1> quads,
    std::pair<const char*, CellDataArrayTypes>... cell_fields)
{
    write_unstructured_grid_header(os, dataset_name, size(quads));
    write_point_coordinates(os, vertices);
    write_cells(os, quads);
    write_cell_data_prelude(os, size(quads));
    (..., write_cell_data(os, cell_fields.first, cell_fields.second));
}

} // namespace vtk
