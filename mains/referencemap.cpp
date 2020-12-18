#include "referencemap.h"

ReferenceMap::ReferenceMap()
{

}

void ReferenceMap::set_velocity_x(CF_2 &velocity_x_in)
{
    velocity_x = &velocity_x_in;
    /// the velocity_x variable in the class is a pointer that points to the address in memory of the velocity_x_in CF_2
}

void ReferenceMap::set_velocity_y(CF_2 &velocity_y_in)
{
    velocity_y = &velocity_y_in;
    /// the velocity_y variable in the class is a pointer that points to the address in memory of the velocity_y_in CF_2
}

void ReferenceMap::set_initial_condition(CF_2 &initial_condition_in)
{
    initial_condition = &initial_condition_in;
}

void ReferenceMap::set_grid(Grid2D &grid_in)
{
    grid = grid_in;
}

void ReferenceMap::set_solution(std::vector<double> &solution_in)
{
    solution = solution_in;
}

void ReferenceMap::set_initial_map_1(std::vector<double> &initial_xi_1_in)
{
    initial_xi_1 = initial_xi_1_in;
}

void ReferenceMap::set_initial_map_2(std::vector<double> &initial_xi_2_in)
{
    initial_xi_2 = initial_xi_2_in;
}

void ReferenceMap::one_Step(double dt)
{
    // do one step for the x component of ref map xi
    std::vector<double> solution_xi_1_n(grid.get_N() * grid.get_M());

    sl_1.set_velocity_x(*velocity_x);
    sl_1.set_velocity_y(*velocity_y);
    sl_1.set_solution(initial_xi_1);
    sl_1.set_grid(grid);

    sl_1.one_Step(dt);
    sl_1.get_solution(solution_xi_1_n);

    // update the ref map x component to be ready for next one step
    initial_xi_1 = solution_xi_1_n;

    // do one step for the y component of ref map xi
    std::vector<double> solution_xi_2_n(grid.get_N() * grid.get_M());

    sl_2.set_velocity_x(*velocity_x);
    sl_2.set_velocity_y(*velocity_y);
    sl_2.set_solution(initial_xi_2);
    sl_2.set_grid(grid);

    sl_2.one_Step(dt);
    sl_2.get_solution(solution_xi_2_n);

    // update the ref map x component to be ready for next one step
    initial_xi_2 = solution_xi_2_n;

    // evaluate initial condition at the reference map
    std::vector<double> solution_n(grid.get_N() * grid.get_M());

#pragma omp parallel for
    for (int n = 0; n < grid.get_N() * grid.get_M(); n++){
        solution_n[n] = (*initial_condition)(solution_xi_1_n[n], solution_xi_2_n[n]);
    }

    set_solution(solution_n);

}

void ReferenceMap::get_solution(std::vector<double> &solution_in)
{
    solution_in = solution;
}















