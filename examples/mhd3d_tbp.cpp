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




#include <iomanip>
#include <sstream>
#include "app_ui.hpp"
#include "core_dimensional.hpp"
#include "core_rational.hpp"
#include "core_sequence.hpp"
#include "core_bqo_tree.hpp"
#include "parallel_dependency_graph.hpp"
#include "parallel_thread_pool.hpp"
#include "scheme_mhd_rules.hpp"




/*
 * Task parallelism should function like a backend to an otherwise functional
 * execution model. There is still a state as before, except that instead of
 * containing the entire solution data, it contains meta-data such as the
 * iteration number, the simulation time, the state of side effects (a.k.a.
 * tasks), the mesh topology, and the run configuration. This state is logically
 * just like the old state, but lazily-evaluated. It *defines* the solution
 * state rather than *contains* it. These "state definitions" map to a list of
 * rules via a lazy version of the advance function. Previously we had advance:
 * state -> state, whereas now we have state_def -> state_def, and we introduce
 * a new function, rules: state_def -> rule_list.
 *
 * This construct is nice in that it separates the flow of states from the
 * workload involved in evaluating them. The logic of profiling, which before
 * was awkwardly stitched into the definition of the simulation state sequence,
 * could be shifted to the runtime, where it probably belongs.
 *
 * The generation of new state definitions requires access to the solution data
 * anytime the mesh topology or time step size depends on the solution data. For
 * this reason we should permit the generator function to read from the graph.
 */




using DependencyGraph = mara::DependencyGraph<mhd_rules::product_identifier_t, mhd_rules::product_t>;
static auto depth = 1;
static auto block_size = 32;




//=============================================================================
void represent_graph_item(std::ostream& os, const DependencyGraph& graph, unsigned row)
{
    auto key = graph.keys().at(row);

    auto shape_string = graph.is_completed(key)
    ? std::visit([] (auto p) { return to_string(p); }, graph.product_at(key))
    : std::string("");

    auto error_string = graph.is_error(key)
    ? std::string("error: ") + graph.error_at(key)
    : std::string("");

    os
    << std::left
    << std::setw(36)
    << std::setfill('.')
    << to_string(key)
    << " status: "
    << to_string(graph.status(key))
    << ' ' << shape_string
    << ' ' << error_string;
}

std::string represent_graph_item(const DependencyGraph& graph, unsigned row)
{
    std::stringstream ss;
    represent_graph_item(ss, graph, row);
    return ss.str();
}




//=============================================================================
auto basic_primitive(mhd_scheme_v2::position_t p, mhd::magnetic_field_vector_t b)
{
    return mhd::primitive(1.0, {}, 1.0, b);
}

auto abc_vector_potential(mhd_scheme_v2::position_t p)
{
    auto k = 2.0 * M_PI / dimensional::unit_length(1.0);
    auto [A, B, C] = std::tuple(1.0, 1.0, 1.0);
    auto [x, y, z] = as_tuple(p);
    auto b0 = 1.0;
    auto ax = A * std::sin(k * z) + C * std::cos(k * y);
    auto ay = B * std::sin(k * x) + A * std::cos(k * z);
    auto az = C * std::sin(k * y) + B * std::cos(k * x);
    return b0 * mhd::vector_potential_t{ax, ay, az};
}




//=============================================================================
struct state_t
{
    rational::number_t                                    iteration;
    dimensional::unit_time                                time;
    dimensional::unit_time                                time_step;
    bsp_tree::shared_tree_t<bqo_tree::tree_index_t<3>, 8> mesh_topology;
};




//=============================================================================
state_t initial_state()
{
    return {
        0,
        0.0,
        0.0,
        bqo_tree::uniform_octree(1),
    };
}

state_t advance(state_t state, dimensional::unit_time dt)
{
    return {
        state.iteration + 1,
        state.time + dt,
        state.time_step = dt,
        state.mesh_topology,
    };
}




//=============================================================================
auto rules(state_t state)
{
    if (long(state.iteration) == 0)
    {
        return flat_map(seq::adapt(state.mesh_topology), [] (auto index)
        {
            return seq::from(mhd_rules::initial_condition(
                index,
                block_size,
                basic_primitive,
                abc_vector_potential));
        }) | seq::to_dynamic();
    }
    else
    {
        return flat_map(seq::adapt(state.mesh_topology), [iteration=state.iteration, dt=state.time_step] (auto index)
        {
            return seq::from(
                mhd_rules::recover_primitive   (index, iteration - 1),
                mhd_rules::extend_primitive    (index, iteration - 1, block_size),
                mhd_rules::extend_magnetic     (index, iteration - 1, block_size),
                mhd_rules::electromotive_force (index, iteration - 1),
                mhd_rules::godunov_data        (index, iteration - 1),
                mhd_rules::update_conserved    (index, iteration, dt),
                mhd_rules::update_magnetic     (index, iteration, dt)
            );
        }) | seq::to_dynamic();
    }
}

auto tasks(state_t state)
{
    return seq::just(std::pair(mhd_rules::product_identifier_t{}, std::function<void(mhd_rules::product_t)>{}));
}




//=============================================================================
class Runtime
{
public:
    using SideEffect    = std::function<void(mhd_rules::product_t)>;
    using SideEffectMap = std::map<mhd_rules::product_identifier_t, SideEffect>;


    //=========================================================================
    Runtime()
    {
        ui_state.runtime_size          = [this] ()               { return graph.size(); };
        ui_state.runtime_item          = [this] (unsigned row)   { return represent_graph_item(graph, row); };
        ui_state.side_effects_size     = [this] ()               { return side_effects.size(); };
        ui_state.side_effects_item     = [this] (unsigned row)   { return to_string(side_effect_at(row)); };
        ui_state.log_view_size         = [this] ()               { return log_vector.size(); };
        ui_state.log_view_item         = [this] (unsigned row)   { return log_vector.at(row); };
        ui_state.concurrent_task_count = [this] ()               { return thread_pool.job_count(); };
        ui::draw(ui_state);
    }

    template<typename RuleIterable, typename TaskIterable>
    void submit(RuleIterable rules, TaskIterable tasks)
    {
        for (auto rule : rules)
        {
            graph.insert(rule);
        }
        for (auto task : tasks)
        {
            side_effects.insert(task);
        }
    }

    void evaluate(std::function<bool(mhd_rules::product_identifier_t)> allow_collection)
    {
        while (graph.count_unevaluated())
        {
            if (ui::fulfill(ui::action::evaluation_step))
            {
                for (const auto& key : graph.eligible_rules())
                {
                    graph.evaluate_rule(key, thread_pool);
                }

                for (const auto& [key, product] : graph.poll(std::chrono::milliseconds(0)))
                {
                    if (side_effects.count(key))
                    {
                        side_effects.at(key)(product);
                        side_effects.erase(key);
                    }
                }

                graph.collect_garbage(allow_collection);
            }

            if (ui::check(ui::action::reset_simulation))
            {
                graph.clear();
                thread_pool.reset(thread_count);
                return;
            }

            if (ui::check(ui::action::quit))
            {
                return;
            }

            if (auto event = ui::poll(std::chrono::milliseconds(20)); event.has_value())
            {
                ui_state = ui::handle_event(ui_state, event.value());
            }
            ui::draw(ui_state);
        }
    }

    const DependencyGraph& get_graph() const
    {
        return graph;
    }

private:

    //=========================================================================
    mhd_rules::product_identifier_t side_effect_at(std::size_t row) const
    {
        for (const auto& s : side_effects)
            if (! row--)
                return s.first;
        throw std::out_of_range("Runtime::side_effect_at");
    }

    //=========================================================================
    unsigned                 thread_count = 1;
    mara::ThreadPool         thread_pool;
    DependencyGraph          graph;
    SideEffectMap            side_effects;
    ui::session_t            termbox;
    ui::state_t              ui_state;
    std::stringstream        log_stream;
    std::vector<std::string> log_vector;
};




//=============================================================================
dimensional::unit_time compute_dt(state_t state, const DependencyGraph& graph)
{
    if (long(state.iteration) > 0)
    {
        // for (auto index : seq::adapt(state.mesh_topology))
        // {
        //     auto key = mhd_rules::key(state.iteration - 1).block(index);
        //     auto uc = std::get<mhd_scheme_v2::cell_conserved_density_t>    (graph.product_at(key.field(mhd_rules::data_field::cell_conserved_density)));
        //     auto bf = std::get<mhd_scheme_v2::face_magnetic_flux_density_t>(graph.product_at(key.field(mhd_rules::data_field::face_magnetic_flux_density)));
        //     mhd_scheme_v2::primitive_array(uc, bf);
        // }
        return dimensional::unit_time(0.001);
    }
    return dimensional::unit_time(0.0);
};




auto maintaining(rational::number_t iteration)
{
    return [iteration] (mhd_rules::product_identifier_t key)
    {
        return ! (
            std::get<0>(key) == iteration && (
            std::get<2>(key) == mhd_rules::data_field::cell_conserved_density ||
            std::get<2>(key) == mhd_rules::data_field::face_magnetic_flux_density));
    };
}




//=============================================================================
int main()
{
    auto runtime = Runtime();
    auto batch_count = 3;

    do
    {
        auto state = initial_state();

        runtime.submit(rules(state), tasks(state));
        runtime.evaluate(maintaining(state.iteration));

        while (true)
        {
            auto dt = compute_dt(state, runtime.get_graph());

            for (int i = 0; i < batch_count; ++i)
            {
                state = advance(state, dt);
                runtime.submit(rules(state), tasks(state));
            }
            runtime.evaluate(maintaining(state.iteration));

            if (ui::fulfill(ui::action::quit))
                return 0;

            if (ui::check(ui::action::reset_simulation))
                break;
        }
    } while (ui::fulfill(ui::action::reset_simulation));

    return 0;
}
