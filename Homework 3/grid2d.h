#ifndef GRID2D_H
#define GRID2D_H

#include <iostream>
#include <vector>

class Grid2D
{
private:
    int N,M;
    double xmin,xmax,ymin,ymax;
    double dx,dy;
public:
    Grid2D();
    Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    double get_dx();
    double get_dy();
    double get_xmin();
    double get_ymin();
    double get_xmax();
    double get_ymax();
    int i_from_n(int n);
    int j_from_n(int n);
    double x_from_n(int n);
    double y_from_n(int n);
    int n_from_ij(int i, int j);

    int get_M();
    int get_N();

    void initialize_VTK_file(std::string file_name);
    void print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name );

    double dx_forward(std::vector<double> &function,int n);
    double dx_backward(std::vector<double> &function,int n);

    double dy_forward(std::vector<double> &function,int n);
    double dy_backward(std::vector<double> &function,int n);

    double sec_dx(std::vector<double> &function,int n);
    double sec_dy(std::vector<double> &function,int n);


    int signum(double val);

};
#endif // GRID2D_H

