/*
#ifndef MATH_TOOLS_H
#define MATH_TOOLS_H

#endif // MATH_TOOLS_H
*/
#include "../../lib/h/grid2d.h"
#include <vector>
#include <cmath>

double bilinear_interpolation(Grid2D &grid, const std::vector<double> &func, const double x, const double y);
double minmod(double f1, double f2);
double abs_val(double f);
double signum(double f, double tol);
double eno_quadratic_interpolation(Grid2D &grid, std::vector<double> &func, const double x, const double y);
