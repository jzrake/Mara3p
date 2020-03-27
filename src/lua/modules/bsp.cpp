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




#define SOL_PRINT_ERRORS 0
#define SOL_ALL_SAFETIES_ON 1
#include "lua/lua.hpp"
#include "lua/sol/sol.hpp"
#include "core_bsp_tree.hpp"
#include "core_bqo_tree.hpp"




//=============================================================================
sol::table open_bsp_lib(sol::this_state s)
{
    using quad_tree_t = bsp::shared_tree<sol::object, 4>;
    using oct_tree_t  = bsp::shared_tree<sol::object, 8>;

    using quad_index_t = bsp::tree_index_t<2>;

    auto lua    = sol::state_view(s);
    auto module = lua.create_table();
    auto quad_tree = module.new_usertype<quad_tree_t>("quad_tree", "", sol::no_constructor);
    auto oct_tree  = module.new_usertype<oct_tree_t> ("oct_tree",  "", sol::no_constructor);

    auto quad_index = module.new_usertype<quad_index_t>("quad_index", sol::call_constructor, sol::default_constructor);

    quad_index["coordinates"] = sol::overload([] (quad_index_t ind) { return ind.coordinates; }, [] (quad_index_t ind, unsigned i, unsigned j) { return bsp::with_coordinates(ind, {i, j}); });
    quad_index["level"]       = sol::overload([] (quad_index_t ind) { return ind.level;       }, [] (quad_index_t ind, unsigned d) { return bsp::with_level(ind, d); });
    quad_index["format"]      = [] (quad_index_t ind) { return bsp::format_tree_index(ind); };
    quad_index[sol::meta_function::to_string] = [] (quad_index_t ind) { return std::string("<bsp::quad_tree_index ") + bsp::format_tree_index(ind) + ">"; };

    quad_tree["index"]   = quad_index;
    quad_tree["just"]    = [] (sol::object value) { return bsp::just<4>(value); };
    quad_tree["uniform"] = [s] (unsigned d) { return bsp::uniform_quadtree(d) | bsp::maps([s] (auto ind) { return sol::make_object(s, ind); }); };

    quad_tree["size"]    = [] (quad_tree_t tree) { return size(tree); };
    quad_tree["start"]   = [] (quad_tree_t tree) { return start(tree); };
    quad_tree["next"]    = [] (quad_tree_t tree, quad_index_t ind) { return next(tree, ind); };
    quad_tree["obtain"]  = [] (quad_tree_t tree, quad_index_t ind) { return obtain(tree, ind); };
    quad_tree["values"]  = [] (quad_tree_t tree) { return std::tuple([] (quad_tree_t tree, quad_index_t ind) { return next(tree, ind); }, tree, start(tree)); };

    quad_tree[sol::meta_function::to_string] = [] (quad_tree_t) { return std::string("<bsp::quad_tree>"); };

    return module;
}
