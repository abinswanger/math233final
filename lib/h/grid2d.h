#ifndef GRID2D_H
#define GRID2D_H

#include <vector>
#include <iostream>

class Grid2D
{
private:
    int N;
    int M;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double dx;
    double dy;
public:
    Grid2D();
    Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    double get_dx();
    double get_dy();
    double get_xmin();
    double get_xmax();
    double get_ymin();
    double get_ymax();
    int get_N();
    int get_M();
    int i_from_n(int n);
    int j_from_n(int n);
    double x_from_n(int n);
    double y_from_n(int n);
    int n_from_ij(int i, int j);
    int number_of_nodes();
    void initialize_VTK_file(std::string file_name);
    void print_VTK_Format(std::vector<double> &F, std::string data_name, std::string file_name);

    // derivatives
    double dx_forward(std::vector<double> &function, int n);
    double dx_backward(std::vector<double> &function, int n);
    double dy_forward(std::vector<double> &function, int n);
    double dy_backward(std::vector<double> &function, int n);
    double dxx_forward(std::vector<double> &function, int n);
    double dxx_backward(std::vector<double> &function, int n);
    double dyy_forward(std::vector<double> &function, int n);
    double dyy_backward(std::vector<double> &function, int n);
    double dxx_central(std::vector<double> &function, int n);
    double dyy_central(std::vector<double> &function, int n);
};

#endif // GRID2D_H
