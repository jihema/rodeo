#include <iostream>
#include "rod_data.h"
#include "rodmath.h"

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
    rod_data.time_step();
    std::cout << rod_data << '\n';
    rod_data.back_step();
    std::cout << rod_data << '\n';

    return 0;
}
