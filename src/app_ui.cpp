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
static bool dummy_session          = false;
static bool wants_quit             = false;
static bool wants_evaluation_step  = false;
static bool wants_reset_simulation = false;
static bool step_evaluations_continuously = false;




//=============================================================================
bool ui::fulfill(action action)
{
    auto check_and_reset = [] (bool& value, bool continuous)
    {
        if (value)
        {
            if (! continuous)
            {
                value = false;
            }
            return true;
        }
        return false;
    };

    switch (action)
    {
        case action::quit            : return check_and_reset(wants_quit, false);
        case action::evaluation_step : return check_and_reset(wants_evaluation_step, step_evaluations_continuously);
        case action::reset_simulation: return check_and_reset(wants_reset_simulation, false);
    }
    return false;
}




//=============================================================================
ui::session_t::session_t(bool is_dummy)
{
    dummy_session = is_dummy;

    if (! dummy_session)
    {
        tb_init();
        tb_select_input_mode(TB_INPUT_ESC | TB_INPUT_MOUSE);
        tb_select_output_mode(TB_OUTPUT_NORMAL);        
    }
    else
    {
        wants_evaluation_step = true;
        step_evaluations_continuously = true;
    }
}

ui::session_t::~session_t()
{
    if (! dummy_session)
    {
        tb_shutdown();        
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
        if (c == '\n')
            c = ' ';
        tb_change_cell(x, y, c, fg, bg);
        x++;
    }
    while (x < min_right)
    {
        tb_change_cell(x, y, ' ', fg, bg);
        x++;        
    }
}

static void draw_usage_tips()
{
    draw_text(right_panel_divider_position() + 2, 6, "play / pause .... p",      TB_BLUE, step_evaluations_continuously ? TB_YELLOW : TB_DEFAULT);
    draw_text(right_panel_divider_position() + 2, 7, "step ............ space",  TB_BLUE, wants_evaluation_step  ? TB_YELLOW : TB_DEFAULT);
    draw_text(right_panel_divider_position() + 2, 8, "reset ........... r",      TB_BLUE, wants_reset_simulation ? TB_YELLOW : TB_DEFAULT);
    draw_text(right_panel_divider_position() + 2, 9, "exit ............ ctrl+q", TB_BLUE, wants_quit             ? TB_YELLOW : TB_DEFAULT);
}

static void draw_navigation_tabs(const ui::state_t& state)
{
    auto tab_names = std::array{"Runtime", "Side Effects", "Log View"};
    auto column = 2;

    for (int i = 0; i < 3; ++i)
    {
        auto bg = state.selected_tab == i ? state.focused_component == ui::component_type::tab_bar
        ? TB_CYAN : TB_BLUE
        : TB_DEFAULT;

        draw_text(column, 2, tab_names.at(i), TB_WHITE | TB_BOLD, bg);
        column += std::strlen(tab_names.at(i)) + 4;
    }
}

static void draw_table(const ui::state_t& state)
{
    std::function<unsigned()>            size_func = nullptr;
    std::function<std::string(unsigned)> item_func = nullptr;
    int starting_row = 0;
    int selected_row = 0;

    switch (state.selected_tab)
    {
        case 0:
            size_func = state.runtime_size;
            item_func = state.runtime_item;
            starting_row = state.starting_runtime_row;
            selected_row = state.selected_runtime_row;
            break;
        case 1:
            size_func = state.side_effects_size;
            item_func = state.side_effects_item;
            starting_row = state.starting_side_effects_row;
            selected_row = state.selected_side_effects_row;
            break;
        case 2:
            size_func = state.log_view_size;
            item_func = state.log_view_item;
            starting_row = state.starting_log_view_row;
            selected_row = state.selected_log_view_row;
            break;
    }

    if (size_func == nullptr) size_func = [] () { return 0; };
    if (item_func == nullptr) item_func = [] (unsigned row) { return "NO DATA"; };

    bool has_focus = state.focused_component == ui::component_type::table;

    for (unsigned i = 0; i < std::min(num_visible_table_rows(), size_func() - starting_row); ++i)
    {
        auto row = starting_row + i;
        auto bg = selected_row == row ? (has_focus ? TB_CYAN : TB_BLUE) : TB_DEFAULT;
        auto fg = item_func(row).find("error") == std::string::npos ? TB_GREEN : TB_RED;

        draw_text(3, 6 + i, item_func(row), fg, bg, right_panel_divider_position() - 1);
    }
}




//=============================================================================
void ui::draw(const state_t& state)
{
    if (dummy_session)
    {
        return;
    }

    int w = tb_width();
    int h = tb_height();

    tb_clear();

    draw_box(0, 0, w, h, TB_CYAN);
    draw_horizontal_line(1, 4, w - 2, TB_CYAN);
    draw_horizontal_line(right_panel_divider_position() + 1, h / 2 + 1, right_panel_width() - 2, TB_CYAN);
    draw_vertical_line(right_panel_divider_position(), 5, h - 6, TB_CYAN);
    draw_navigation_tabs(state);
    draw_usage_tips();
    draw_text(right_panel_divider_position() + 2, h / 2 + 2, "Job count: " + std::to_string(state.concurrent_task_count()), TB_YELLOW);

    draw_table(state);
    tb_present();
}




//=============================================================================
bool ui::is_quit(tb_event ev)
{
    return ev.key == TB_KEY_CTRL_Q;
}

bool ui::is_dummy_session()
{
    return dummy_session;
}




//=============================================================================
ui::state_t move_focus_forward(ui::state_t state)
{
    switch (state.focused_component)
    {
         case ui::component_type::tab_bar:       state.focused_component = ui::component_type::table; break;
         case ui::component_type::table:         state.focused_component = ui::component_type::detail_view; break;
         case ui::component_type::detail_view:   state.focused_component = ui::component_type::options_panel; break;
         case ui::component_type::options_panel: state.focused_component = ui::component_type::tab_bar; break;
    }
    return state;
}

static ui::state_t scroll_table(ui::state_t state, int difference)
{
    auto scroll = [difference] (unsigned& starting, unsigned size)
    {
        if (int(starting) + difference + num_visible_table_rows() <= size &&
            int(starting) + difference >= 0)
        {
            starting += difference;
        }
    };
    switch (state.selected_tab)
    {
        case 0: scroll(state.starting_runtime_row,      state.runtime_size()); break;
        case 1: scroll(state.starting_side_effects_row, state.side_effects_size()); break;
        case 2: scroll(state.starting_log_view_row,     state.log_view_size()); break;
    }
    return state;
}

static ui::state_t select_prev_table_item(ui::state_t state)
{
    auto dec = [] (unsigned& selected, unsigned& starting)
    {
        if (selected > 0)
        {
            selected -= 1;

            if (selected < starting)
            {
                starting -= 1;
            }
        }
    };
    switch (state.selected_tab)
    {
        case 0: dec(state.selected_runtime_row,      state.starting_runtime_row); break;
        case 1: dec(state.selected_side_effects_row, state.starting_side_effects_row); break;
        case 2: dec(state.selected_log_view_row,     state.starting_log_view_row); break;
    }
    return state;
}

static ui::state_t select_next_table_item(ui::state_t state)
{
    auto inc = [&state] (unsigned& selected, unsigned& starting)
    {
        unsigned size = 0;

        switch (state.selected_tab)
        {
            case 0: size = state.runtime_size      ? state.runtime_size()      : 0; break;
            case 1: size = state.side_effects_size ? state.side_effects_size() : 0; break;
            case 2: size = state.log_view_size     ? state.log_view_size()     : 0; break;
        }

        if (selected < size - 1)
        {
            selected += 1;

            if (selected >= starting + num_visible_table_rows())
            {
                starting += 1;
            }
        }
    };
    switch (state.selected_tab)
    {
        case 0: inc(state.selected_runtime_row,      state.starting_runtime_row); break;
        case 1: inc(state.selected_side_effects_row, state.starting_side_effects_row); break;
        case 2: inc(state.selected_log_view_row,     state.starting_log_view_row); break;
    }
    return state;
}

static ui::state_t select_table_at_position(ui::state_t state, int y)
{
    auto sel = [y] (unsigned& selected, unsigned starting)
    {
        selected = y - 6 + starting;
    };
    switch (state.selected_tab)
    {
        case 0: sel(state.selected_runtime_row,      state.starting_runtime_row); break;
        case 1: sel(state.selected_side_effects_row, state.starting_side_effects_row); break;
        case 2: sel(state.selected_log_view_row,     state.starting_log_view_row); break;
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
    if (ev.ch == 'p')     step_evaluations_continuously = ! step_evaluations_continuously;

    if (step_evaluations_continuously)
    {
        wants_evaluation_step = true;
    }
    return state;
}

ui::state_t handle_mouse(tb_event ev, ui::state_t state)
{
    switch (ev.key)
    {
        case TB_KEY_MOUSE_LEFT:       return select_table_at_position(state, ev.y);
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
