#include "../../lib/h/sl_method.h"

SL_method::SL_method()
{

}

void SL_method::set_velocity_x(CF_2 &velocity_x_in)
{
    velocity_x = &velocity_x_in;
    /// the velocity_x variable in the class is a pointer that points to the address in memory of the velocity_x_in CF_2
}

void SL_method::set_velocity_y(CF_2 &velocity_y_in)
{
    velocity_y = &velocity_y_in;
    /// the velocity_y variable in the class is a pointer that points to the address in memory of the velocity_y_in CF_2
}

void SL_method::set_grid(Grid2D &grid_in)
{
    grid = grid_in;
}

void SL_method::set_solution(std::vector<double> &solution_in)
{
    solution = solution_in;
}

std::vector<double> SL_method::trajectory_interpolation(int n, double dt)
{
    // RK2
    double xstar = grid.x_from_n(n) - 0.5*dt*(*velocity_x)(grid.x_from_n(n),grid.y_from_n(n));
    double ystar = grid.y_from_n(n) - 0.5*dt*(*velocity_y)(grid.x_from_n(n),grid.y_from_n(n));

    double xd = grid.x_from_n(n) - dt*(*velocity_x)(xstar,ystar);
    double yd = grid.y_from_n(n) - dt*(*velocity_y)(xstar,ystar);

    std::vector<double> dep_pt(2);
    dep_pt[0] = xd;
    dep_pt[1] = yd;

    return dep_pt;
}

void SL_method::one_Step(double dt)
{
#pragma omp parallel for
    for(int n = 0; n < grid.number_of_nodes(); n++){
        std::vector<double> dep_pt = trajectory_interpolation(n,dt);
        solution[n] = eno_quadratic_interpolation(grid, solution, dep_pt[0], dep_pt[1]);
    }
}

void SL_method::get_solution(std::vector<double> &solution_in)
{
    solution_in = solution;
}

