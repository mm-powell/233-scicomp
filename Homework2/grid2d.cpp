#include "grid2d.h"
#include <vector>

Grid2D::Grid2D()
{
}
Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
{
    N = N_;
    M = M_;
    xmin = xmin_;
    xmax = xmax_;
    ymin = ymin_;
    ymax = ymax_;
    dx = (xmax - xmin) / (double) (N - 1);
    dy = (ymax - ymin) / (double) (M - 1);
}

double Grid2D::get_dx()
{
    return dx;
}

double Grid2D::get_dy()
{
    return dy;
}

int Grid2D::get_M()
{
    return M;
}

int Grid2D::get_N()
{
    return N;
}

int Grid2D::i_from_n(int n)
{
    // n = i+ j * N
    return  n % N;
}
int Grid2D::j_from_n(int n)
{
    // n = i+ j * N
    return  n / N;
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




// Number 4
double Grid2D::dx_forward(std::vector<double> &function,int n)
{
    if ( int(function.size()) != (N*M) ) { std::cout<<"Error!"<<std::endl; return 1.; }
    // else if ( x_from_n(n) == xmax ) {return (function[n] - function[n-1])/dx;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( y_from_n(n) == ymin ) {return 0.0;}
    else if ( y_from_n(n) == ymax ) {return 0.0;}
    else {return (function[n+1] - function[n])/dx;}
}

double Grid2D::dx_backward(std::vector<double> &function,int n)
{

    if ( int(function.size()) != (N*M) ) { std::cout<<"Error!"<<std::endl; return 1.; }
    // else if ( x_from_n(n) == xmin ) {return (function[n+1] - function[n])/dx;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( y_from_n(n) == ymin ) {return 0.0;}
    else if ( y_from_n(n) == ymax ) {return 0.0;}
    else {return (function[n] - function[n-1])/dx;}
}

// Number 5
double Grid2D::dy_forward(std::vector<double> &function,int n)
{
    if ( int(function.size()) != (N*M) ) { std::cout<<"Error!"<<std::endl; return 1.; }
    // else if ( y_from_n(n) == ymax ) {return (function[n] - function[n-1])/dy;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( y_from_n(n) == ymin ) {return 0.0;}
    else if ( y_from_n(n) == ymax ) {return 0.0;}
    else { return (function[n+N] - function[n])/dy; }
}

double Grid2D::dy_backward(std::vector<double> &function,int n)
{
    if ( int(function.size()) != (N*M) ) { std::cout<<"Error!"<<std::endl; return 1.; }
    // else if ( y_from_n(n) == ymin ) {return (function[n+1] - function[n])/dy;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( x_from_n(n) == xmax ) {return 0.0;}
    else if ( y_from_n(n) == ymin ) {return 0.0;}
    else if ( y_from_n(n) == ymax ) {return 0.0;}
    else {return (function[n] - function[n-N])/dy;}
}


// Number 6
double Grid2D::sec_dx(std::vector<double> &function,int n) //use first order derivative to approx second order derivative
{
    // throw an exception
    if ( int(function.size()) != (N*M) )
    {
        std::cout<<"Error!"<<std::endl; return 1.;
    }

    if ( x_from_n(n) == xmin )
    {
         return 0.;
        //return ((function[n] - (2.0*function[n+1]) + function[n+2])/(dx*dx));
    }
    else if ( x_from_n(n) == xmax )
    {
         return 0.;
        //return ((function[n] - (2.0*function[n-1]) + function[n-2])/(dx*dx));
    }
    //else {return ((function[n+1] - (2*function[n]) + function[n-1])/(dx*dx));}
//    else if ( x_from_n(n) == xmin ) { return 0.; }
//    else if ( x_from_n(n) == xmax ) { return 0.; }
    //
    else {return ((function[n+1] - (2.0*function[n]) + function[n-1])/(dx*dx));}
}

double Grid2D::sec_dy(std::vector<double> &function,int n)
{
    if ( int(function.size()) != (N*M) )
    {
        std::cout<<"Error!"<<std::endl; return 1.;
    }
//    else if ( y_from_n(n) == ymin ) { return ((function[n] - 2*function[n+1] + function[n+2])/(dy*dy)); }
//    else if ( y_from_n(n) == ymax ) { return ((function[n-2] - 2*function[n-1] + function[n])/(dy*dy)); }

    if ( y_from_n(n) == ymin )
    {
        return 0.;
        //return ((function[n] - (2.0*function[n+N]) + function[n+(2*N)])/(dy*dy));
    }

    else if ( y_from_n(n) == ymax || y_from_n(n-N) == ymax  )
    {
        return 0.;
        //return ((function[n-(2*N)] - (2.0*function[n-N]) + function[n])/(dy*dy));
    }

//    else if ( y_from_n(n) == ymin ) { return 0.; }
//    else if ( y_from_n(n) == ymax ) { return 0.; }
    //else {return ((function[n+1] - (2*function[n]) + function[n-1])/(dy*dy));}

    else {

//        std::cout<<n+N << std::endl;

        return ((function[n+N] - (2.0*function[n]) + function[n-N])/(dy*dy));
    }
}
