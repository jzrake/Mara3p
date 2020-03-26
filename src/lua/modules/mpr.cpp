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
#include <sstream>
#include "lua/lua.hpp"
#include "lua/sol/sol.hpp"
#include "parallel_computable.hpp"




//=============================================================================
sol::table open_mpr_lib(sol::this_state s)
{
    using lua_computable = mpr::computable<sol::object>;

    auto lua        = sol::state_view(s);
    auto module     = lua.create_table();
    auto computable = module.new_usertype<lua_computable>("computable", "", sol::no_constructor);

    computable["has_value"] = &lua_computable::has_value;
    computable["value"]     = &lua_computable::value;
    computable["map"]       = [] (lua_computable c, sol::function f) { return mpr::map(c, f); };
    computable["compute"]   = [] (lua_computable c) { return mpr::compute(c); };
    computable["as_dot"]    = [] (lua_computable c) { auto stream = std::stringstream(); mpr::print_graph(stream, c); return stream.str(); };
    computable["name"]      = sol::overload([] (lua_computable c) { return c.name(); }, [] (lua_computable c, std::string name) { return c.name(name); });
    computable["immediate"] = sol::overload([] (lua_computable c) { return c.immediate(); }, [] (lua_computable c, bool imm) { return c.immediate(imm); });

    module["just"] = mpr::just<sol::object>;
    module["from"] = [] (sol::function f) { return mpr::from([f] () -> sol::object { return f(); }); };
    module["zip"]  = [] (sol::this_state s, sol::variadic_args args)
    {
        auto nodes = mpr::node_set_t();

        for (auto a : args)
        {
            nodes.insert(a.as<lua_computable>().node());
        }

        return lua_computable([args, s] ()
        {
            auto lua = sol::state_view(s);
            auto result = lua.create_table();
            auto n = int(0);

            for (auto a : args)
            {
                result[n++] = a.as<lua_computable>().value();
            }
            return result;
        }, nodes).immediate(true);
    };

    return module;
}
