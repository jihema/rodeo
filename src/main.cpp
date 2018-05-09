#include <iostream>
#include "rod_data.h"

int main()
{
    static constexpr size_t num_vertices = 15;

    rodeo::RodData rod_data(num_vertices);

    std::cout << rod_data << '\n';
    rod_data.time_step();
    std::cout << rod_data << '\n';
    rod_data.back_step();
    std::cout << rod_data << '\n';

    return 0;
}
