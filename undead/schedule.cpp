#include <iostream>
#include "sequence.hpp"
#include "numeric_tuple.hpp"
#include "numeric_array.hpp"




//=============================================================================
struct task_t
{
    std::string name;
    unsigned long count = 0;
    double next_time = 0.0;
};

auto jump(task_t task, double time, double cadence)
{
    if (time < task.next_time)
    {
        task.count += 1;
        task.next_time += cadence;
    }
    return task;
}

auto jump_task(double time, std::function<double(std::string)> cadence)
{
    return [=] (task_t task)
    {
        return jump(task, time, cadence(task.name));
    };
}




//=============================================================================
struct state_t
{
    int iteration = 0;
    double time = 0.0;
    numeric::array_t<task_t, 3> tasks;
};




//=============================================================================
state_t initial_state()
{
    return state_t{
        0,
        0.0,
        {
            task_t{"task1"},
            task_t{"task2"},
            task_t{"task3"}
        },
    };
}

state_t next_state(state_t s)
{
    auto cadence = [] (std::string name)
    {
        if (name == "task1") return 0.10;
        if (name == "task2") return 0.25;
        if (name == "task3") return 0.91;
        throw std::invalid_argument("no such task");
    };

    return {
        s.iteration + 1,
        s.time + 0.102301,
        map(s.tasks, jump_task(s.time, cadence))
    };
}

bool should_continue(std::pair<state_t, state_t> state)
{
    return state.first.time < 2.0;
}

auto side_effects(std::pair<state_t, state_t> states)
{
    return [s=states.second] ()
    {
        std::printf("[%04d]: t=%1.6lf\n", s.iteration, s.time);
    };
}




//=============================================================================
int main()
{
    auto simulation = seq::generate(initial_state(), next_state)
    | seq::window()
    | seq::take_while(should_continue)
    | seq::map(side_effects);

    for (auto side_effect : simulation)
    {
        side_effect();
    }
    return 0;
}

