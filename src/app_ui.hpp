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
#include <chrono>
#include <functional>
#include <optional>
#include <string>
#include <vector>
#include <termbox.h>




//=============================================================================
namespace ui
{


//=============================================================================
enum class component_type
{
    tab_bar,
    content_table,
    detail_view,
    options_panel,
};

enum class action
{
    quit,
    evaluation_step,
    reset_simulation,
};




//=============================================================================
struct state_t
{
    component_type focused_component;

    unsigned selected_tab = 0;
    unsigned selected_table_row = 0;
    unsigned starting_table_row = 0;

    std::function<unsigned()>            concurrent_task_count;
    std::function<unsigned()>            content_table_size;
    std::function<std::string(unsigned)> content_table_item;
};

struct session_t
{
    session_t(bool is_dummy_session=false);
    ~session_t();
};




//=============================================================================
bool is_quit(tb_event ev);
bool is_dummy_session();
void draw(const state_t& state);
bool fulfill(action action);

ui::state_t handle_event(const state_t& state, tb_event ev);

template<class Rep, class Period>
std::optional<tb_event> poll(std::chrono::duration<Rep, Period> timeout)
{
    if (! is_dummy_session())
    {
        auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(timeout);
        tb_event ev;

        if (tb_peek_event(&ev, dt.count()) > 0)
        {
            return ev;
        }
    }
    return {};
}

} // namespace ui
