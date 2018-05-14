/*
 * explicit_euler.cpp
 */

#include "explicit_euler.h"
#include "../rods/rod_dynamic.h"
#include "../forces/force_base.h"

namespace rodeo
{

void ExplicitEulerSolver::compute(force::ForceBase * force,
        ExplicitEulerRodState& current_state)
{
    if (!force)
    {
        return;
    }

    force->accumulate_energy_force(&current_state);
}

void ExplicitEulerSolver::time_step(ExplicitEulerRodState& future_state,
        ExplicitEulerRodState const& current_state, double const h)
{
    assert(!current_state.is_dirty());

    try
    {

        // Copy the (pointers to) fixed parameters.
        static_cast<RodBase&>(future_state) = current_state;

        VecXd const acceleration = current_state.get_inverse_mass()
                * current_state.get_force();

        // Explicit step in velocity.
        future_state.set_velocity(
                current_state.get_velocity() + h * acceleration);

        // Explicit step in dofs.
        future_state.set_dofs(&current_state.get_reference_frames(),
                current_state.get_dofs() + h * current_state.get_velocity()
                        + 0.5 * sqr(h) * acceleration);

        // Record the time change.
        future_state.set_time(current_state.get_time() + h);
    }
    catch (...) // Something went wrong, just copy current to future (including time).
                // It is thus possible that in the future times will be inconsistent.
    {
        std::cerr << "Warning: in " << __PRETTY_FUNCTION__
                << " exception caught, evolution skipped.";
        future_state = current_state;
    }
}

}

