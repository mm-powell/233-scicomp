#include "math_tools.h"
#include <math.h>
#include <cmath>

double minmod( double vel1, double vel2 )
{
    if ( vel1*vel2 < 0.0) {return 0.0;}
    else if ( std::abs(vel1) < std::abs(vel2) ) { return vel1; }
    else { return vel2; }
}



double bilinear_interpolation(Grid2D & grid, std::vector<double> & func, const double x, const double y)
{
    int i_min = (int) floor( (x - grid.get_xmin()) / grid.get_dx() );
        i_min = std::max(0, i_min);
    int i_max = (int) ceil ( (x - grid.get_xmin()) / grid.get_dx() );
        i_max = std::min(grid.get_N() - 1, i_max);


    if(i_min == i_max) {
        if(i_min==0){
            i_max = i_min+1;
        } else{
            i_min = i_max-1;
        }
    }


    int j_min = (int) floor( (y - grid.get_ymin()) / grid.get_dy() );
        j_min = std::max(0, j_min);
    int j_max = (int) ceil ( (y - grid.get_ymin()) / grid.get_dy() );
        j_max = std::min(grid.get_M() - 1, j_max);

    if(j_min == j_max) {
        if(j_min==0){
            j_max = j_min+1;
        } else{
            j_min = j_max-1;
        }
    }


    int corn_LB = grid.n_from_ij(i_min, j_min);
    int corn_LT = grid.n_from_ij(i_min, j_max);
    int corn_RB = grid.n_from_ij(i_max, j_min);
    int corn_RT = grid.n_from_ij(i_max, j_max);

    double dx = grid.get_dx();
    double dy = grid.get_dy();

    double x_min = grid.x_from_n(corn_LB);
    double y_min = grid.y_from_n(corn_LB);

    double x_max = grid.x_from_n(corn_RT);
    double y_max = grid.y_from_n(corn_RT);


    return ( 1.0 / ( dx * dy ) ) * ( func[corn_LB] * ( x_max - x ) * ( y_max - y ) +
                                     func[corn_LT] * ( x_max - x ) * ( y - y_min ) +
                                     func[corn_RB] * ( x - x_min ) * ( y_max - y ) +
                                     func[corn_RT] * ( x - x_min ) * ( y - y_min ) );
}









double ENO_interpolation(Grid2D & grid, std::vector<double> & func, const double x, const double y)
{
    int i_min = (int) floor( (x - grid.get_xmin()) / grid.get_dx() );
        i_min = std::max(0, i_min);
    int i_max = (int) ceil ( (x - grid.get_xmin()) / grid.get_dx() );
        i_max = std::min(grid.get_N() - 1, i_max);


    if(i_min == i_max) {
        if(i_min==0){
            i_max = i_min+1;
        } else{
            i_min = i_max-1;
        }
    }

    int j_min = (int) floor( (y - grid.get_ymin()) / grid.get_dy() );
        j_min = std::max(0, j_min);
    int j_max = (int) ceil ( (y - grid.get_ymin()) / grid.get_dy() );
        j_max = std::min(grid.get_M() - 1, j_max);

    if(j_min == j_max) {
        if(j_min==0){
            j_max = j_min+1;
        } else{
            j_min = j_max-1;
        }
    }


    int corn_00 = grid.n_from_ij(i_min, j_min);
    int corn_01 = grid.n_from_ij(i_min, j_max);
    int corn_10 = grid.n_from_ij(i_max, j_min);
    int corn_11 = grid.n_from_ij(i_max, j_max);

    double dx = grid.get_dx();
    double dy = grid.get_dy();

    double x_min = grid.x_from_n(corn_00);
    double y_min = grid.y_from_n(corn_00);

    double x_max = grid.x_from_n(corn_11);
    double y_max = grid.y_from_n(corn_11);


    double dx2_00 = grid.sec_dx(func, corn_00);
    double dx2_01 = grid.sec_dx(func, corn_01);
    double dx2_10 = grid.sec_dx(func, corn_10);
    double dx2_11 = grid.sec_dx(func, corn_11);

    double xxmm_a = minmod(dx2_00, dx2_01);
    double xxmm_b = minmod(dx2_10, dx2_11);
    double xx_minmod = minmod(xxmm_a, xxmm_b);

    double dy2_00 = grid.sec_dy(func, corn_00);
    double dy2_01 = grid.sec_dy(func, corn_01);
    double dy2_10 = grid.sec_dy(func, corn_10);
    double dy2_11 = grid.sec_dy(func, corn_11);

    double yymm_a = minmod(dy2_00, dy2_01);
    double yymm_b = minmod(dy2_10, dy2_11);
    double yy_minmod = minmod(yymm_a, yymm_b);


    return (  (( 1.0 / ( dx * dy ) ) * ( func[corn_00] * ( x_max - x ) * ( y_max - y ) +
                                         func[corn_01] * ( x_max - x ) * ( y - y_min ) +
                                         func[corn_10] * ( x - x_min ) * ( y_max - y ) +
                                         func[corn_11] * ( x - x_min ) * ( y - y_min ) ))
            -  ( (0.5)*( (x - x_min)*(x_max - x) )*( xx_minmod ) )
            -  ( (0.5)*( (y - y_min)*(y_max - y) )*( yy_minmod ) )  );
}

