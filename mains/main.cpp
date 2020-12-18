#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

#include "../../lib/h/sl_method.h"
#include "../../lib/h/math_tools.h"
#include "../../lib/h/grid2d.h"
#include <../../lib/h/cf_2.h>
#include "referencemap.h"

using namespace std;

double t_sol = 0.;

class velocity_X :public CF_2
{
public:
    double operator()(double x, double y) const{
        return -y + 0.*x;
    }
}velocity_x;

class velocity_Y :public CF_2
{
public:
    double operator()(double x, double y) const{
        return x + 0.*y;
    }
}velocity_y;

class initial_Condition :public CF_2
{
public:
    double operator()(double x, double y) const{
        if( (y > -0.2) && (y < (-2.*x + 1.2) ) && (y < (2.*x - 0.8) ) && (y < 0.2) )
        {
            return 1.;
        }
        return 0.;
    }
}initial_condition;

class initial_Map_1 :public CF_2
{
public:
    double operator()(double x, double y) const{
        return x + 0.*y;
    }
}initial_map_1;

class initial_Map_2 :public CF_2
{
public:
    double operator()(double x, double y) const{
        return y + 0.*x;
    }
}initial_map_2;

class true_Solution :public CF_2
{
public:
    double operator()(double x, double y) const{
        if( (y*cos(t_sol) - x*sin(t_sol) > -0.2) && (y*cos(t_sol) - x*sin(t_sol) < (-2.*(x*cos(t_sol) + y*sin(t_sol)) + 1.2) ) && (y*cos(t_sol) - x*sin(t_sol) < (2.*(x*cos(t_sol) + y*sin(t_sol)) - 0.8) ) && (y*cos(t_sol) - x*sin(t_sol) < 0.2) )
        {
            return 1.;
        }
        return 0.;
    }
}true_solution;


int main()
{
    int N = 128;
    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;

    Grid2D grid(N,N,xmin,xmax,ymin,ymax);

    // vector of initial condition values defined at the nodes of the grid
    vector<double> initial_condition_n(N*N);

    // vector of reference map initial conditions defined at the nodes of the grid
    vector<double> initial_map_1_n(N*N);
    vector<double> initial_map_2_n(N*N);

    // vector of true solution
    vector<double> true_solution_n(N*N);

#pragma omp parallel for
    for (int n=0; n<N*N; n++){
        initial_condition_n[n] = initial_condition(grid.x_from_n(n),grid.y_from_n(n));
        initial_map_1_n[n] = initial_map_1(grid.x_from_n(n),grid.y_from_n(n));
        initial_map_2_n[n] = initial_map_2(grid.x_from_n(n),grid.y_from_n(n));
        true_solution_n[n] = true_solution(grid.x_from_n(n),grid.y_from_n(n));
    }

    // set semi lagrangian solver parameters
    SL_method sl;
    sl.set_velocity_x(velocity_x);
    sl.set_velocity_y(velocity_y);
    sl.set_solution(initial_condition_n);
    sl.set_grid(grid);

    // set reference map class parameters
    ReferenceMap rf;
    rf.set_velocity_x(velocity_x);
    rf.set_velocity_y(velocity_y);
    rf.set_initial_condition(initial_condition);
    rf.set_initial_map_1(initial_map_1_n);
    rf.set_initial_map_2(initial_map_2_n);
    rf.set_grid(grid);

    // define time step and number of associated steps
    double dtdxratio = 1. / 5.;

    double dt = dtdxratio*grid.get_dx();
    double tf = 2*M_PI;

    int nmax = floor(tf/dt);

    // add initial profiles to vtk file
    char name[250];
    sprintf(name,"/home/adam/Documents/math233vtk/refmap2/grid2D_N=%d_M=%d_%d.vtk",N,N,0);
    grid.initialize_VTK_file(name);
    grid.print_VTK_Format(initial_condition_n,"solution_sl",name);
    grid.print_VTK_Format(initial_condition_n,"solution_rf",name);
    grid.print_VTK_Format(true_solution_n,"true_solution",name);

    // define vectors to populate solution at each time step
    vector<double> solution_sl_n(N*N);
    vector<double> solution_rf_n(N*N);

    double t = dt;
    t_sol = t;

    for(int n = 1; n < nmax; n++)
    {
        for (int k=0; k<N*N; k++){
            true_solution_n[k] = true_solution(grid.x_from_n(k),grid.y_from_n(k));
        }

        sl.one_Step(dt);
        rf.one_Step(dt);
        sprintf(name,"/home/adam/Documents/math233vtk/refmap2/grid2D_N=%d_M=%d_%d.vtk",N,N,n);
        grid.initialize_VTK_file(name);
        sl.get_solution(solution_sl_n);
        rf.get_solution(solution_rf_n);
        grid.print_VTK_Format(solution_sl_n,"solution_sl",name);
        grid.print_VTK_Format(solution_rf_n,"solution_rf",name);
        grid.print_VTK_Format(true_solution_n,"true_solution",name);
        cout<<"printing in vtk file at "<<n<<"th time step"<<endl;

        t += dt;
        t_sol = t;
    }

    // adapt the time step and get solution at final time
    double dt_final = tf - t;
    sl.one_Step(dt_final);
    rf.one_Step(dt_final);

    sprintf(name,"/home/adam/Documents/math233vtk/refmap2/grid2D_N=%d_M=%d_%d.vtk",N,N,nmax+1);
    grid.initialize_VTK_file(name);
    sl.get_solution(solution_sl_n);
    rf.get_solution(solution_rf_n);
    grid.print_VTK_Format(solution_sl_n,"solution_sl",name);
    grid.print_VTK_Format(solution_rf_n,"solution_rf",name);

    double inferror = 0.;
    double inferrortemp = 0.;

    t_sol = tf;

    for (int k=0; k<N*N; k++){
        true_solution_n[k] = true_solution(grid.x_from_n(k),grid.y_from_n(k));

        inferrortemp = abs(true_solution_n[k] - solution_rf_n[k]);
        if(inferrortemp > inferror){
            inferror = inferrortemp;
        }
    }

    cout<<"inf error "<<inferror<<endl;

    cout << "Hello World!" << endl;
    return 0;
}















