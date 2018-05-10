/*
 * explicit_euler.cpp
 */

#include "explicit_euler.h"
#include "rod_state.h"

namespace rodeo
{

void ExplicitEulerSolver::time_step(ExplicitEulerRodState& future_state, ExplicitEulerRodState const& current_state,
        double const h)
{
    assert(!current_state.is_dirty());

    // Copy the dofs and frames from current state to future.
    // FIXME: technically we don't need to copy material frames...
    static_cast<StaticRodState<0>&>(future_state) = current_state;

    future_state.set_dofs(current_state.get_dofs() + h * current_state.get_velocity(), false);

    future_state.set_velocity(
            current_state.get_velocity() + h * current_state.get_inverse_mass() * current_state.get_force());
}

}

