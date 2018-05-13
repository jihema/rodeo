#include "forces/gravitation.h"
#include "forces/stretching.h"
#include "forces/composite.h"
#include "rods/rod_data.h"
#include "rodmath.h"
#include "solvers/explicit_euler.h"

#include <iostream>

int main()
{
    using namespace rodeo;

    static constexpr size_t num_vertices = 15;

    VecXd rest_dofs(4 * num_vertices - 1);
    for (size_t i = 0; i < num_vertices; ++i)
    {
        rest_dofs[4 * i] = i;
    }

    RodData rod_data(rest_dofs);

    force::Gravitation gravitation;
    force::Stretching stretching(1.);

    force::Composite all_force;
    all_force.add_component(&gravitation);
    all_force.add_component(&stretching);
    rod_data.set_force(&all_force);

    for (int i = 0; i < 15; ++i)
    {
        rod_data.forward_step(0.01);
    }

    return 0;
}
