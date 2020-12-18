#include "../../lib/h/math_tools.h"
#include "../../lib/h/grid2d.h"

#include <vector>
#include <math.h>
#include <cmath>

#include <iostream>


double bilinear_interpolation(Grid2D &grid, const std::vector<double> & func, const double x, const double y)
{
    // extract indices of four corners of cell that (x,y) lies in
    int imin = (int) floor((x - grid.get_xmin())/grid.get_dx());
    imin = std::max(0,imin); // left of domain
    imin = std::min(imin,grid.get_N()-2); // right of domain
    int imax = (int) ceil((x - grid.get_xmin())/grid.get_dx());
    //int imax = imin++;
    imax = std::min(grid.get_N()-1,imax); // right of domain
    imax = std::max(1,imax); // left of domain

    if(imin == imax){
        if(imin == 0){
            imax = imin + 1;
        }else{
            imin = imax - 1;
        }
    }

    int jmin = (int) floor((y - grid.get_ymin())/grid.get_dy());
    jmin = std::max(0,jmin); // below domain
    jmin = std::min(jmin,grid.get_M()-2); // above domain
    int jmax = (int) ceil((y - grid.get_ymin())/grid.get_dy());
    //int jmax = jmin++;
    jmax = std::min(grid.get_M()-1,jmax); // above domain
    jmax = std::max(1,jmax); // below domain

    if(jmin == jmax){
        if(jmin == 0){
            jmax = jmin + 1;
        }else{
            jmin = jmax - 1;
        }
    }

    // find n values from above indices, grid has bottom left corner 00 and top right 11
    int n_00 = grid.n_from_ij(imin,jmin);
    int n_10 = grid.n_from_ij(imax,jmin);
    int n_01 = grid.n_from_ij(imin,jmax);
    int n_11 = grid.n_from_ij(imax,jmax);

    // find corresponding x and y values
    double xi = grid.x_from_n(n_00);
    double yj = grid.y_from_n(n_00);
    double xip1 = grid.x_from_n(n_11);
    double yjp1 = grid.y_from_n(n_11);

    // bilinear interpolation
    double interpval = func[n_00]*(xip1 - x)*(yjp1 - y);
    interpval += func[n_01]*(xip1 - x)*(y - yj);
    interpval += func[n_10]*(x - xi)*(yjp1 - y);
    interpval += func[n_11]*(x - xi)*(y - yj);
    interpval *= 1./(grid.get_dx()*grid.get_dy());

    return interpval;
}

double minmod(double f1, double f2)
{
    if(f1*f2 > 0){
        return std::min(abs(f1),abs(f2));
    }
    return 0.;
}

double signum(double f, double tol)
{

    if(f > 0.){
        if(f < tol){
            return 0.;
        }else{
            return 1.;
        }
    }

    if(f < 0.){
        if(f > -tol){
            return 0.;
        }else{
            return -1.;
        }
    }

    return 0.;
}

double eno_quadratic_interpolation(Grid2D &grid, std::vector<double> &func, const double x, const double y)
{
    // extract indices of four corners of cell that (x,y) lies in
    int imin = (int) floor((x - grid.get_xmin())/grid.get_dx());
    imin = std::max(0,imin); // left of domain
    imin = std::min(imin,grid.get_N()-2); // right of domain
    int imax = (int) ceil((x - grid.get_xmin())/grid.get_dx());
    //int imax = imin++;
    imax = std::min(grid.get_N()-1,imax); // right of domain
    imax = std::max(1,imax); // left of domain

    if(imin == imax){
        if(imin == 0){
            imax = imin + 1;
        }else{
            imin = imax - 1;
        }
    }

    int jmin = (int) floor((y - grid.get_ymin())/grid.get_dy());
    jmin = std::max(0,jmin); // below domain
    jmin = std::min(jmin,grid.get_M()-2); // above domain
    int jmax = (int) ceil((y - grid.get_ymin())/grid.get_dy());
    //int jmax = jmin++;
    jmax = std::min(grid.get_M()-1,jmax); // above domain
    jmax = std::max(1,jmax); // below domain

    if(jmin == jmax){
        if(jmin == 0){
            jmax = jmin + 1;
        }else{
            jmin = jmax - 1;
        }
    }

    int boundary_checker = 0;

    // check if on wall for eno second derivatives
    if(imin == 0){
        // left wall
        boundary_checker = 1;
    }

    if(jmin == 0){
        // bottom wall
        boundary_checker = 1;
    }

    if(imax == grid.get_N()-1){
        // right wall
        boundary_checker = 1;
    }

    if(jmax == grid.get_M()-1){
        // top wall
        boundary_checker = 1;
    }

    // find n values from above indices, grid has bottom left corner 00 and top right 11
    int n_00 = grid.n_from_ij(imin,jmin);
    int n_10 = grid.n_from_ij(imax,jmin);
    int n_01 = grid.n_from_ij(imin,jmax);
    int n_11 = grid.n_from_ij(imax,jmax);

    // find corresponding x and y values
    double xi = grid.x_from_n(n_00);
    double yj = grid.y_from_n(n_00);
    double xip1 = grid.x_from_n(n_11);
    double yjp1 = grid.y_from_n(n_11);

    // bilinear interpolation
    double interpval = func[n_00]*(xip1 - x)*(yjp1 - y);
    interpval += func[n_01]*(xip1 - x)*(y - yj);
    interpval += func[n_10]*(x - xi)*(yjp1 - y);
    interpval += func[n_11]*(x - xi)*(y - yj);
    interpval *= 1./(grid.get_dx()*grid.get_dy());

    // add quadratic (eno inspired) term if not on boundary
    if(boundary_checker == 0){
        interpval -= 0.5*(x - xi)*(xip1 - x)*minmod(grid.dxx_central(func,n_00),minmod(grid.dxx_central(func,n_01),minmod(grid.dxx_central(func,n_10),grid.dxx_central(func,n_11))));
        interpval -= 0.5*(y - yj)*(yjp1 - y)*minmod(grid.dyy_central(func,n_00),minmod(grid.dyy_central(func,n_01),minmod(grid.dyy_central(func,n_10),grid.dyy_central(func,n_11))));
    }

    return interpval;

}












