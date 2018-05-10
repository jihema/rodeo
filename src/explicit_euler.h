/*
 * explicit_euler.h
 */

#pragma once


namespace rodeo
{

class ExplicitEulerRodState;

class ExplicitEulerSolver
{
public:

    static void time_step(ExplicitEulerRodState& future_state, ExplicitEulerRodState const& current_state,
            double const h);

};

}

