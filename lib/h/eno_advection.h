#ifndef ENO_ADVECTION_H
#define ENO_ADVECTION_H

#include "../../lib/h/grid2d.h"
#include <vector>
#include <../../lib/h/cf_2.h>
#include <cmath>
#include<omp.h>

class ENO_Advection
{
public:
    ENO_Advection();
    void set_velocity_x(CF_2 &velocity_x_in);
    void set_velocity_y(CF_2 &velocity_y_in);
    void set_solution(std::vector<double> &solution_in);
    void set_grid(Grid2D &grid_in);
    void one_Step(double dt);
    void get_solution(std::vector<double> &solution_in);
    void one_Step_central(double dt);

    void print_v(); //for debugging

private:
    Grid2D eno_grid;
    std::vector<double> solution;
    CF_2 *velocity_x;
    CF_2 *velocity_y;
    double velocity_x_value(double x, double y);
    double velocity_y_value(double x, double y);
    double minmod(double f1, double f2);

};

#endif // ENO_ADVECTION_H
