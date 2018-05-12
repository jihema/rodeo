#include "forces/gravitation.h"
#include "forces/stretching.h"
#include "forces/composite.h"
#include "rods/rod_data.h"
#include "rodmath.h"
#include "solvers/explicit_euler.h"

#include <iostream>

int main()
{
    static constexpr size_t num_vertices = 15;

    rodeo::VecXd rest_dofs(4*num_vertices-1);
    for (size_t i = 0; i < num_vertices;++i)
    {
        rest_dofs[4*i]=i;
    }

    rodeo::RodData rod_data(rest_dofs);

    std::cout << rod_data << '\n';
    rod_data.time_step(0.1);
    std::cout << rod_data << '\n';
    rod_data.time_step(0.1);
    std::cout << rod_data << '\n';

    return 0;
}
