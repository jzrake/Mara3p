/**
 ==============================================================================
 Copyright 2020, Jonathan Zrake

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




#include "mara.hpp"
#if MARA_COMPILE_LOCALLY_ISOTHERMAL
#include "app_hdf5_dimensional.hpp"
#include "app_hdf5_ndarray.hpp"
#include "app_hdf5_ndarray_dimensional.hpp"
#include "app_hdf5_numeric_array.hpp"
#include "core_bqo_tree.hpp"
#include "core_util.hpp"
#include "module_locally_isothermal.hpp"
#endif




#if MARA_COMPILE_LOCALLY_ISOTHERMAL // <---------------------------------------


//=============================================================================
void h5::read(const Group& group, std::string name,  modules::locally_isothermal::conserved_tree_t& conserved)
{
    auto mesh = bsp::just<4>(mesh::block_index_t<2>());

    for (auto block_name : group)
    {
        auto block = mesh::read_block_index<2>(block_name);
        mesh = insert(mesh, block, block);
    }

    conserved = mesh | bsp::maps([&group] (auto block)
    {
        return mpr::from([block, &group] ()
        {
            return read<modules::locally_isothermal::conserved_array_t>(group, mesh::format_block_index(block));
        });
    });
}




//=============================================================================
void h5::write(const Group& group, std::string name, const modules::locally_isothermal::conserved_tree_t& conserved)
{
    auto ugroup = group.require_group(name);

    sink(indexify(conserved), util::apply_to([&ugroup] (auto index, auto value)
    {
        if (value.has_value())
        {
            write(ugroup, mesh::format_block_index(index), value.value());
        }
    }));
}

#endif // MARA_COMPILE_LOCALLY_ISOTHERMAL
