/*
 * explicit_euler.cpp
 */

#include "explicit_euler.h"
#include "../rods/rod_dynamic.h"

namespace rodeo
{

void ExplicitEulerSolver::time_step(ExplicitEulerRodState& future_state,
        ExplicitEulerRodState const& current_state, double const h)
{
    assert(!current_state.is_dirty());

    // Copy the (pointers to) fixed parameters.
    static_cast<RodBase&>(future_state) = current_state;

    VecXd const acceleration = current_state.get_inverse_mass()
            * current_state.get_force();

    // Record the time change.
    future_state.set_time(current_state.get_time() + h);

    // Explicit step in velocity.
    future_state.set_velocity(current_state.get_velocity() + h * acceleration);

    // Explicit step in dofs.
    future_state.set_dofs(&current_state.get_reference_frames(),
            current_state.get_dofs() + h * current_state.get_velocity()
                    + 0.5 * sqr(h) * acceleration);
}

}

