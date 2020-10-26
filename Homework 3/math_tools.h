#include <grid2d.h>
#include <vector>

double bilinear_interpolation(Grid2D & grid, std::vector<double> & func, const double x, const double y);
double ENO_interpolation(Grid2D & grid, std::vector<double> & func, const double x, const double y);
double eno_interp(Grid2D & grid, std::vector<double> & func, const double x, const double y);
double minmod( double vel1, double vel2 );

double reinit(Grid2D & grid, std::vector<double> & func, const double x, const double y);
