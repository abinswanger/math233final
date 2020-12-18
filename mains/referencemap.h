#ifndef REFERENCEMAP_H
#define REFERENCEMAP_H

#include "../../lib/h/grid2d.h"
#include <vector>
#include <../../lib/h/cf_2.h>
#include "../../lib/h/math_tools.h"
#include "../../lib/h/sl_method.h"
#include "omp.h"

class ReferenceMap
{
private:
    Grid2D grid;
    std::vector<double> solution;
    std::vector<double> initial_xi_1;
    std::vector<double> initial_xi_2;

    CF_2 *velocity_x;
    CF_2 *velocity_y;
    CF_2 *initial_condition;

    SL_method sl_1;
    SL_method sl_2;

    void set_solution(std::vector<double> &solution_in);

public:
    ReferenceMap();
    void set_velocity_x(CF_2 &velocity_x_in);
    void set_velocity_y(CF_2 &velocity_y_in);
    void set_initial_condition(CF_2 &initial_condition_in);
    void set_initial_map_1(std::vector<double> &initial_xi_1_in);
    void set_initial_map_2(std::vector<double> &initial_xi_2_in);
    void set_grid(Grid2D &grid_in);
    void one_Step(double dt);
    void get_solution(std::vector<double> &solution_in);
};

#endif // REFERENCEMAP_H
