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




#include <array>
#include "app_ui.hpp"




//=============================================================================
static bool wants_quit             = false;
static bool wants_evaluation_step  = false;
static bool wants_reset_simulation = false;




//=============================================================================
bool ui::fulfill(action action)
{
    auto check_and_reset = [] (bool& value)
    {
        if (value)
        {
            value = false;
            return true;
        }
        return false;
    };

    switch (action)
    {
        case action::quit            : return check_and_reset(wants_quit);
        case action::evaluation_step : return check_and_reset(wants_evaluation_step);
        case action::reset_simulation: return check_and_reset(wants_reset_simulation);
    }
    return false;
}




//=============================================================================
ui::session_t::session_t()
{
    tb_init();
    tb_select_input_mode(TB_INPUT_ESC | TB_INPUT_MOUSE);
}

ui::session_t::~session_t()
{
    tb_shutdown();
}




//=============================================================================
static void draw_box(int x, int y, int w, int h, uint16_t fg=TB_WHITE, uint16_t bg=TB_DEFAULT)
{
    int x0 = x;
    int x1 = x + w - 1;
    int y0 = y;
    int y1 = y + h - 1;

    tb_change_cell(x0, y0, 0x250C, fg, bg);
    tb_change_cell(x1, y0, 0x2510, fg, bg);
    tb_change_cell(x0, y1, 0x2514, fg, bg);
    tb_change_cell(x1, y1, 0x2518, fg, bg);

    for (int i = x0 + 1; i < x1; ++i)
    {
        tb_change_cell(i, y0, 0x2500, fg, bg);
        tb_change_cell(i, y1, 0x2500, fg, bg);
    }

    for (int i = y0 + 1; i < y1; ++i)
    {
        tb_change_cell(x0, i, 0x2502, fg, bg);
        tb_change_cell(x1, i, 0x2502, fg, bg);
    }
}

static void draw_horizontal_line(int x, int y, int w, uint16_t fg=TB_WHITE, uint16_t bg=TB_DEFAULT)
{
    for (int i = x; i < x + w; ++i)
    {
        tb_change_cell(i, y, 0x2500, fg, bg);
    }
}

static void draw_vertical_line(int x, int y, int h, uint16_t fg=TB_WHITE, uint16_t bg=TB_DEFAULT)
{
    for (int j = y; j < y + h; ++j)
    {
        tb_change_cell(x, j, 0x2502, fg, bg);
    }
}

static void draw_text(int x, int y, const std::string& text, uint16_t fg=TB_WHITE, uint16_t bg=TB_DEFAULT, unsigned min_right=0)
{
    for (auto c : text)
    {
        tb_change_cell(x, y, c, fg, bg);
        x++;
    }
    while (x < min_right)
    {
        tb_change_cell(x, y, ' ', fg, bg);
        x++;        
    }
}




//=============================================================================
static unsigned num_visible_table_rows()
{
    return tb_height() - 8;
}

static unsigned right_panel_width()
{
    return 40;
}

static unsigned right_panel_divider_position()
{
    return tb_width() - right_panel_width();
}




//=============================================================================
void ui::draw(const state_t& state)
{
    int w = tb_width();
    int h = tb_height();

    tb_clear();
    draw_box(0, 0, w, h, TB_CYAN);
    draw_horizontal_line(1, 4, w - 2, TB_CYAN);
    draw_horizontal_line(right_panel_divider_position() + 1, h / 2 + 1, right_panel_width() - 2, TB_CYAN);
    draw_vertical_line(right_panel_divider_position(), 5, h - 6, TB_CYAN);




    auto tab_names = std::array{"Execution Graph", "Performance", "File Browser"};
    auto column = 3;

    for (int i = 0; i < 3; ++i)
    {
        auto bg = state.selected_tab == i ? state.focused_component == component_type::tab_bar
        ? TB_CYAN : TB_BLUE
        : TB_DEFAULT;

        draw_text(column, 2, tab_names.at(i), TB_WHITE | TB_BOLD, bg);
        column += std::strlen(tab_names.at(i)) + 4;
    }


    draw_text(right_panel_divider_position() + 2, 6, "step ...... space",  TB_BLUE, wants_evaluation_step  ? TB_YELLOW : TB_DEFAULT);
    draw_text(right_panel_divider_position() + 2, 7, "reset ..... r",      TB_BLUE, wants_reset_simulation ? TB_YELLOW : TB_DEFAULT);
    draw_text(right_panel_divider_position() + 2, 8, "exit ...... ctrl+q", TB_BLUE, wants_quit             ? TB_YELLOW : TB_DEFAULT);


    draw_text(right_panel_divider_position() + 2, h / 2 + 2, "Job count: " + std::to_string(state.concurrent_task_count));



    if (state.selected_tab == 0)
    {
        for (unsigned i = 0; i < num_visible_table_rows(); ++i)
        {
            auto row = state.starting_table_row + i;

            if (row >= state.content_table_items.size())
            {
                break;
            }

            auto bg = state.selected_table_row == row
            ? state.focused_component == component_type::content_table ? TB_CYAN : TB_BLUE
            : TB_DEFAULT;

            auto fg = state.content_table_items.at(row).find("error") == std::string::npos ? TB_GREEN : TB_RED;

            draw_text(3, 6 + i, state.content_table_items.at(row), fg, bg, 120);
        }
    }
    tb_present();
}




//=============================================================================
std::optional<tb_event> ui::peek(int timeout)
{
    tb_event ev;

    if (tb_peek_event(&ev, timeout) > 0)
    {
        return ev;
    }
    return {};
}




//=============================================================================
bool ui::is_quit(tb_event ev)
{
    return ev.key == TB_KEY_CTRL_Q;
}




//=============================================================================
ui::state_t move_focus_forward(ui::state_t state)
{
    switch (state.focused_component)
    {
         case ui::component_type::tab_bar:       state.focused_component = ui::component_type::content_table; break;
         case ui::component_type::content_table: state.focused_component = ui::component_type::detail_view; break;
         case ui::component_type::detail_view:   state.focused_component = ui::component_type::options_panel; break;
         case ui::component_type::options_panel: state.focused_component = ui::component_type::tab_bar; break;
    }
    return state;
}

static ui::state_t scroll_table(ui::state_t state, int difference)
{
    if (int(state.starting_table_row) + difference + num_visible_table_rows() <= state.content_table_items.size() &&
        int(state.starting_table_row) + difference >= 0)
    {
        state.starting_table_row += difference;
        return state;
    }
    return state;
}

static ui::state_t select_prev_table_item(ui::state_t state)
{
    if (state.selected_table_row > 0)
    {
        state.selected_table_row -= 1;

        if (state.selected_table_row < state.starting_table_row)
        {
            state.starting_table_row -= 1;
        }
    }
    return state;
}

static ui::state_t select_next_table_item(ui::state_t state)
{
    if (state.selected_table_row < state.content_table_items.size() - 1)
    {
        state.selected_table_row += 1;

        if (state.selected_table_row >= state.starting_table_row + num_visible_table_rows())
        {
            state.starting_table_row += 1;
        }
    }
    return state;
}




//=============================================================================
static ui::state_t handle_key(tb_event ev, ui::state_t state)
{
    switch (ev.key)
    {
    case TB_KEY_ARROW_UP     : return select_prev_table_item(state);
    case TB_KEY_ARROW_DOWN   : return select_next_table_item(state);
    case TB_KEY_ARROW_LEFT   : if (state.selected_tab > 0) state.selected_tab -= 1; return state;
    case TB_KEY_ARROW_RIGHT  : if (state.selected_tab < 2) state.selected_tab += 1; return state;
    case TB_KEY_TAB          : return move_focus_forward(state);
    case TB_KEY_SPACE        : wants_evaluation_step = true; return state;
    }

    if (ui::is_quit(ev))  wants_quit = true;
    if (ev.ch == 'r')     wants_reset_simulation = true;

    return state;
}

ui::state_t handle_mouse(tb_event ev, ui::state_t state)
{
    switch (ev.key)
    {
        case TB_KEY_MOUSE_WHEEL_UP:   return scroll_table(state, -1);
        case TB_KEY_MOUSE_WHEEL_DOWN: return scroll_table(state, +1);
    }
    return state;
}




//=============================================================================
ui::state_t ui::handle_event(const ui::state_t& state, tb_event ev)
{
    switch (ev.type)
    {
    case TB_EVENT_KEY:    return handle_key(ev, state);
    case TB_EVENT_MOUSE:  return handle_mouse(ev, state);
    case TB_EVENT_RESIZE: return state;
    }
    return state;
}
