#include "../../lib/h/eno_advection.h"
#include <iostream>
#include<omp.h>
using namespace std;

ENO_Advection::ENO_Advection()
{

}

void ENO_Advection::set_velocity_x(CF_2 &velocity_x_in)
{
    velocity_x = &velocity_x_in;
    /// the velocity_x variable in the class is a pointer that points to the address in memory of the velocity_x_in CF_2
}

void ENO_Advection::set_velocity_y(CF_2 &velocity_y_in)
{
    velocity_y = &velocity_y_in;
    /// the velocity_y variable in the class is a pointer that points to the address in memory of the velocity_y_in CF_2
}

double ENO_Advection::velocity_x_value(double x, double y)
{
    return velocity_x->operator ()(x,y);
}

double ENO_Advection::velocity_y_value(double x, double y)
{
    return velocity_y->operator ()(x,y);
}

void ENO_Advection::set_grid(Grid2D &grid_in)
{
    eno_grid = grid_in;
}

void ENO_Advection::set_solution(std::vector<double> &solution_in)
{
    solution = solution_in;
}

void ENO_Advection::get_solution(std::vector<double> &solution_in)
{
    solution_in = solution;
}

double ENO_Advection::minmod(double f1, double f2)
{
    if(f1*f2 > 0){
        return min(abs(f1),abs(f2));
    }
    return 0.;
}

void ENO_Advection::one_Step(double dt)
{
#pragma omp parallel for
    for(int n = 0; n < eno_grid.number_of_nodes(); n++){
        int i = eno_grid.i_from_n(n);
        int j = eno_grid.j_from_n(n);

        int im = i-1;
        int jm = j-1;

        int ip = i+1;
        int jp = j+1;

        int nmx = eno_grid.n_from_ij(im,j);
        int npx = eno_grid.n_from_ij(ip,j);

        int nmy = eno_grid.n_from_ij(i,jm);
        int npy = eno_grid.n_from_ij(i,jp);

        double x_n = eno_grid.x_from_n(n);
        double y_n = eno_grid.y_from_n(n);

        double u_n = velocity_x_value(x_n,y_n);
        double v_n = velocity_y_value(x_n,y_n);

        int boundary_checker = 0;

        // check if on the boundary first


        if((i == 0) && (j == 0)){
            // bottom left corner

            /*
            solution[n] += u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,npx);
            solution[n] += v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,npy);*/
            boundary_checker = 1;
        }

        if((j == 0) && (i > 0) && (i < eno_grid.get_N()-1)){
            // bottom wall but not a corner
            /*
            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n) + 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,nmx))
                                                 : eno_grid.dx_forward(solution,n)  - 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,npx)));
            solution[n] += v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,npy);*/
            boundary_checker = 1;
        }

        if((i == eno_grid.get_N()-1) && (j == 0)){
            // bottom right corner
            /*
            solution[n] -= u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,nmx);
            solution[n] += v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,npy);*/
            boundary_checker = 1;
        }

        if((i == 0) && (j > 0) && (j < eno_grid.get_M()-1)){
            // left wall but not a corner
            /*
            solution[n] += u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,npx);
            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n) + 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,nmy))
                                                 : eno_grid.dy_forward(solution,n)  - 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,npy)));*/
            boundary_checker = 1;
        }

        if((i == 0) && (j == eno_grid.get_M()-1)){
            // top left corner
            /*
            solution[n] += u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,npx);
            solution[n] -= v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,nmy);*/
            boundary_checker = 1;
        }

        if((j == eno_grid.get_M()-1) && (i > 0) && (i < eno_grid.get_N()-1)){
            // top wall but not a corner
            /*
            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n) + 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,nmx))
                                                 : eno_grid.dx_forward(solution,n)  - 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,npx)));
            solution[n] -= v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,nmy);*/
            boundary_checker = 1;
        }

        if((i == eno_grid.get_N()-1) && (j == eno_grid.get_M()-1)){
            // top right corner
            /*
            solution[n] -= u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,nmx);
            solution[n] -= v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,nmy);*/
            boundary_checker = 1;
        }

        if((i == eno_grid.get_N()-1) && (j > 0) && (j < eno_grid.get_M()-1)){
            // right wall but not a corner
            /*
            solution[n] -= u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,nmx);
            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n) + 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,nmy))
                                                 : eno_grid.dy_forward(solution,n)  - 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,npy)));*/
            boundary_checker = 1;
        }

        if(boundary_checker == 0){
            // interior points

            /*
            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n)
                                                 : eno_grid.dx_forward(solution,n));

            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n)
                                                 : eno_grid.dy_forward(solution,n));*/

            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n) + 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,nmx))
                                                 : eno_grid.dx_forward(solution,n)  - 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,npx)));

            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n) + 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,nmy))
                                                 : eno_grid.dy_forward(solution,n)  - 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,npy)));
        }
    }
}

void ENO_Advection::one_Step_central(double dt)
{
#pragma omp parallel for
    for(int n = 0; n < eno_grid.number_of_nodes(); n++){
        int i = eno_grid.i_from_n(n);
        int j = eno_grid.j_from_n(n);

        double x_n = eno_grid.x_from_n(n);
        double y_n = eno_grid.y_from_n(n);

        double u_n = velocity_x_value(x_n,y_n);
        double v_n = velocity_y_value(x_n,y_n);

        int boundary_checker = 0;

        // check if on the boundary first


        if((i == 0) && (j == 0)){
            // bottom left corner

            /*
            solution[n] += u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,npx);
            solution[n] += v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,npy);*/
            boundary_checker = 1;
        }

        if((j == 0) && (i > 0) && (i < eno_grid.get_N()-1)){
            // bottom wall but not a corner
            /*
            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n) + 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,nmx))
                                                 : eno_grid.dx_forward(solution,n)  - 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,npx)));
            solution[n] += v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,npy);*/
            boundary_checker = 1;
        }

        if((i == eno_grid.get_N()-1) && (j == 0)){
            // bottom right corner
            /*
            solution[n] -= u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,nmx);
            solution[n] += v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,npy);*/
            boundary_checker = 1;
        }

        if((i == 0) && (j > 0) && (j < eno_grid.get_M()-1)){
            // left wall but not a corner
            /*
            solution[n] += u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,npx);
            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n) + 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,nmy))
                                                 : eno_grid.dy_forward(solution,n)  - 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,npy)));*/
            boundary_checker = 1;
        }

        if((i == 0) && (j == eno_grid.get_M()-1)){
            // top left corner
            /*
            solution[n] += u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,npx);
            solution[n] -= v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,nmy);*/
            boundary_checker = 1;
        }

        if((j == eno_grid.get_M()-1) && (i > 0) && (i < eno_grid.get_N()-1)){
            // top wall but not a corner
            /*
            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n) + 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,nmx))
                                                 : eno_grid.dx_forward(solution,n)  - 0.5*eno_grid.get_dx()*minmod(eno_grid.dxx_central(solution,n),eno_grid.dxx_central(solution,npx)));
            solution[n] -= v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,nmy);*/
            boundary_checker = 1;
        }

        if((i == eno_grid.get_N()-1) && (j == eno_grid.get_M()-1)){
            // top right corner
            /*
            solution[n] -= u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,nmx);
            solution[n] -= v_n * dt * 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,nmy);*/
            boundary_checker = 1;
        }

        if((i == eno_grid.get_N()-1) && (j > 0) && (j < eno_grid.get_M()-1)){
            // right wall but not a corner
            /*
            solution[n] -= u_n * dt * 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,nmx);
            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n) + 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,nmy))
                                                 : eno_grid.dy_forward(solution,n)  - 0.5*eno_grid.get_dy()*minmod(eno_grid.dyy_central(solution,n),eno_grid.dyy_central(solution,npy)));*/
            boundary_checker = 1;
        }

        if(boundary_checker == 0){
            // interior points

            /*
            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n)
                                                 : eno_grid.dx_forward(solution,n));

            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n)
                                                 : eno_grid.dy_forward(solution,n));*/

            solution[n] -= u_n * dt * ( u_n > 0. ? eno_grid.dx_backward(solution,n) + 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,n)
                                                 : eno_grid.dx_forward(solution,n)  - 0.5*eno_grid.get_dx()*eno_grid.dxx_central(solution,n));

            solution[n] -= v_n * dt * ( v_n > 0. ? eno_grid.dy_backward(solution,n) + 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,n)
                                                 : eno_grid.dy_forward(solution,n)  - 0.5*eno_grid.get_dy()*eno_grid.dyy_central(solution,n));
        }
    }
}

void ENO_Advection::print_v()
{
    cout<<velocity_x<<endl;
    cout<<velocity_y<<endl;

    double val = velocity_x_value(1,1);
    double val2 = velocity_y_value(1,1);

    cout<<val<<endl;
    cout<<val2<<endl;

    cout<<eno_grid.get_dx()<<endl;

    //string food = "pizza";
    //string *ptr = &food;
    //cout<<ptr<<endl;
    //cout<<*ptr<<endl;

}
