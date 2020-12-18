#include "../../lib/h/grid2d.h"
#include <iostream>

Grid2D::Grid2D()
{

}

Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_){
    N = N_;
    M = M_;
    xmin = xmin_;
    xmax = xmax_;
    ymin = ymin_;
    ymax = ymax_;

    dx = (xmax-xmin) / (double) (N-1);
    dy = (ymax-ymin) / (double) (M-1);
}

double Grid2D::get_dx()
{
    return dx;
}

double Grid2D::get_dy()
{
    return dy;
}

double Grid2D::get_xmin()
{
    return xmin;
}

double Grid2D::get_xmax()
{
    return xmax;
}

double Grid2D::get_ymin()
{
    return ymin;
}

double Grid2D::get_ymax()
{
    return ymax;
}

int Grid2D::get_N()
{
    return N;
}

int Grid2D::get_M()
{
    return M;
}

int Grid2D::i_from_n(int n){
    // n = i + j*N
    return n % N;
}

int Grid2D::j_from_n(int n){
    // n = i + j*N
    return n / N;
}

int Grid2D::n_from_ij(int i, int j)
{
    return i+ j * N;
}

double Grid2D::x_from_n(int n)
{
    return  xmin + dx*i_from_n(n);
}
double Grid2D::y_from_n(int n)
{
    return  ymin + dy*j_from_n(n);
}

double Grid2D::dx_forward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int ip = i+1;

    int np = n_from_ij(ip,j);

    return (function[np]-function[n])/dx;
}

double Grid2D::dx_backward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int im = i-1;

    int nm = n_from_ij(im,j);

    return (function[n]-function[nm])/dx;
}

double Grid2D::dy_forward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int jp = j+1;

    int np = n_from_ij(i,jp);

    return (function[np]-function[n])/dy;
}

double Grid2D::dy_backward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int jm = j-1;

    int nm = n_from_ij(i,jm);

    return (function[n]-function[nm])/dy;
}

double Grid2D::dxx_backward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int im = i-1;
    int imm = i-2;

    int nm = n_from_ij(im,j);
    int nmm = n_from_ij(imm,j);

    return (function[n] - 2*function[nm] + function[nmm]) / (dx*dx);
}

double Grid2D::dxx_forward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int ip = i+1;
    int ipp = i+2;

    int np = n_from_ij(ip,j);
    int npp = n_from_ij(ipp,j);

    return (function[npp] - 2*function[np] + function[n]) / (dx*dx);
}

double Grid2D::dyy_backward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int jm = j-1;
    int jmm = j-2;

    int nm = n_from_ij(i,jm);
    int nmm = n_from_ij(i,jmm);

    return (function[n] - 2*function[nm] + function[nmm]) / (dy*dy);
}

double Grid2D::dyy_forward(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int jp = j+1;
    int jpp = j+2;

    int np = n_from_ij(i,jp);
    int npp = n_from_ij(i,jpp);

    return (function[npp] - 2*function[np] + function[n]) / (dy*dy);
}

double Grid2D::dxx_central(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int ip = i+1;
    int im = i-1;

    int np = n_from_ij(ip,j);
    int nm = n_from_ij(im,j);

    return (function[np] - 2*function[n] + function[nm]) / (dx*dx);
}

double Grid2D::dyy_central(std::vector<double> &function, int n)
{
    if((int) function.size() != N*M){
        throw "e";
    }

    int i = i_from_n(n);
    int j = j_from_n(n);
    int jp = j+1;
    int jm = j-1;

    int np = n_from_ij(i,jp);
    int nm = n_from_ij(i,jm);

    return (function[np] - 2*function[n] + function[nm]) / (dy*dy);
}

int Grid2D::number_of_nodes()
{
    return N*M;
}

// initialize the .vtk file at the specified address with all the grid information

void Grid2D::initialize_VTK_file(std::string file_name)
{
  int node_of_cell[4];

  FILE *outFile = fopen(file_name.c_str(),"w");

  fprintf(outFile,"# vtk DataFile Version 2.0 \n");
  fprintf(outFile,"Quadtree Mesh \n");
  fprintf(outFile,"ASCII \n");
  fprintf(outFile,"DATASET UNSTRUCTURED_GRID \n");


//% first output the list of nodes
  fprintf(outFile,"POINTS %d double \n",N*M);
  for (int n=0; n<N*M; n++)
    fprintf(outFile,"%e %e %e\n",x_from_n(n), y_from_n(n), 0.0);


  // then output the list of cells. each cell is composed of four nodes
  fprintf(outFile,"CELLS %d %d \n",(N-1)*(M-1),5*(N-1)*(M-1));
  for (int i=0; i<N-1; i++)
      for (int j=0; j<M-1; j++)
      {
          node_of_cell[0] = n_from_ij(i  ,j  );
          node_of_cell[1] = n_from_ij(i+1,j  );
          node_of_cell[2] = n_from_ij(i+1,j+1);
          node_of_cell[3] = n_from_ij(i  ,j+1);

          fprintf(outFile,"%d %d %d %d %d\n",4,node_of_cell[0], node_of_cell[1], node_of_cell[2], node_of_cell[3]);
          }
  //  }
  fprintf(outFile,"CELL_TYPES %d \n",(N-1)*(M-1));
  for (int n=0; n<(N-1)*(M-1); n++)    fprintf(outFile,"%d \n",9);
  fprintf(outFile,"POINT_DATA %d \n",N*M);
  fclose (outFile);
}

// this function write the values of the vector F into the vtk file. before using it, the .vtk file must have been initialized with all the grid infos
void Grid2D::print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name )
{

  FILE *outFile;
  outFile = fopen(file_name.c_str(),"a");
  fprintf(outFile,"SCALARS %s double 1 \n",data_name.c_str());
  fprintf(outFile,"LOOKUP_TABLE default \n");
  for (int n=0; n<N*M; n++) fprintf(outFile,"%e \n",F[n]);
  fclose (outFile);
}
