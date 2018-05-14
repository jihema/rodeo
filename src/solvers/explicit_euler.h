/*
 * explicit_euler.h
 */

#pragma once

namespace rodeo
{

namespace force
{
class ForceBase;
}

class ExplicitEulerRodState;

class ExplicitEulerSolver
{
public:
    static void compute(force::ForceBase* force,
            ExplicitEulerRodState& current_state);

    static void time_step(ExplicitEulerRodState& future_state,
            ExplicitEulerRodState const& current_state, double const h);
};

}

