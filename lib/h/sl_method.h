#ifndef SL_METHOD_H
#define SL_METHOD_H

#include "../../lib/h/grid2d.h"
#include <vector>
#include <../../lib/h/cf_2.h>
#include "../../lib/h/math_tools.h"

class SL_method
{
private:
    Grid2D grid;
    std::vector<double> solution;
    CF_2 *velocity_x;
    CF_2 *velocity_y;
    std::vector<double> trajectory_interpolation(int n, double dt);

public:
    SL_method();
    void set_velocity_x(CF_2 &velocity_x_in);
    void set_velocity_y(CF_2 &velocity_y_in);
    void set_solution(std::vector<double> &solution_in);
    void set_grid(Grid2D &grid_in);
    void one_Step(double dt);
    void get_solution(std::vector<double> &solution_in);

};

#endif // SL_METHOD_H
