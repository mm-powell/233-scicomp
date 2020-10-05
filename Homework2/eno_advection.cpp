#include "eno_advection.h"
#include <grid2d.h>
#include <cmath>
#include <omp.h>
#include <iostream>

ENO_Advection::ENO_Advection()
{

}

ENO_Advection::ENO_Advection(Grid2D grid_, std::vector<double> ans_vec_, velocity_X * velocityx_, velocity_Y * velocityy_)
{
    grid = grid_;
    ans_vec = ans_vec_;
    velocityx = velocityx_;
    velocityy = velocityy_;
}

double ENO_Advection::minmod( double vel1, double vel2 )
{
    if ( vel1*vel2 < 0.0) {return 0.0;}
    else if ( std::abs(vel1) < std::abs(vel2) ) { return vel1; }
    else { return vel2; }
}

std::vector<double> ENO_Advection::get_sol()
{
    return ans_vec;
}

void ENO_Advection::Central_Advection(double dt)
{
    for (int n=0; n < (grid.get_N()*grid.get_M()) ; ++n)
    {
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);

        double dx = grid.get_dx();
        double dy = grid.get_dy();

//        double vx = velocityx->operator ()(x,y);
//        double vy = velocityy->operator ()(x,y);
        double vx = (* velocityx)(x,y);
        double vy = (* velocityy)(x,y);

//        if ( vx > 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( ((vx) * ( grid.dx_backward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) )) +
//                                                                         ((vy) * ( grid.dy_backward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) )) ); }
//        else if ( vx < 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( vx * ( grid.dx_forward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) ) +
//                                                                              vy * ( grid.dy_forward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) ));}
//        else if ( vx < 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( vx * ( grid.dx_forward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) ) +
//                                                                              vy * ( grid.dy_backward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) )); }
//        else if ( vx > 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( vx * ( grid.dx_backward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) ) +
//                                                                              vy * ( grid.dy_forward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) )); }
//        else {ans_vec[n] = ans_vec[n];}

//        if ( vx > 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( ((vx) * ( grid.dx_backward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) )) +
//                                                                         ((vy) * ( grid.dy_backward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) )) ); }
//        else if ( vx < 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( vx * ( grid.dx_forward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) ) +
//                                                                              vy * ( grid.dy_forward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) ));}
//        else if ( vx < 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( vx * ( grid.dx_forward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) ) +
//                                                                              vy * ( grid.dy_backward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) )); }
//        else if ( vx > 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( vx * ( grid.dx_backward(ans_vec, n) + (grid.sec_dx(ans_vec, n)*dx) ) +
//                                                                              vy * ( grid.dy_forward(ans_vec, n) + (grid.sec_dy(ans_vec, n)*dy) )); }
//        else {ans_vec[n] = ans_vec[n];}

        if ( vx >= 0.0 && vy >= 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_backward(ans_vec, n) + ( grid.sec_dx(ans_vec,n) )*0.5*dx))  +
                                                                            (vy * ( grid.dy_backward(ans_vec, n) + (( grid.sec_dy(ans_vec,n) )*0.5*dy))) ) ; }
        else if ( vx < 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_forward(ans_vec, n) + (( grid.sec_dx(ans_vec,n) )*0.5*dx)))  +
                                                                                 (vy * ( grid.dy_forward(ans_vec, n) + (( grid.sec_dy(ans_vec,n) )*0.5*dy))) ); }
        else if ( vx < 0.0 && vy >= 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_forward(ans_vec, n) + (( grid.sec_dx(ans_vec,n) )*0.5*dx)))  +
                                                                                 (vy * ( grid.dy_backward(ans_vec, n) + ( grid.sec_dy(ans_vec,n) )*0.5*dy))) ; }
        else if ( vx >= 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_backward(ans_vec, n) + (( grid.sec_dx(ans_vec,n) )*0.5*dx)))  +
                                                                                 (vy * ( grid.dy_forward(ans_vec, n) + (( grid.sec_dy(ans_vec,n) )*0.5*dy))) ); }


    }
}


void ENO_Advection::one_Step(double dt)
{
    #pragma omp parallel for
    for (int n=0; n < (grid.get_N()*grid.get_M()) ; ++n)
    {
        int N = grid.get_N();
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);

        double dx = grid.get_dx();
        double dy = grid.get_dy();

        double vx = velocityx->operator ()(x,y);
        double vy = velocityy->operator ()(x,y);
        //double vx = (* velocityx)(x,y);
        //double vy = (* velocityy)(x,y);


//        if ( vx > 0.0 ) {
//            dx_ = grid.dx_backward(ans_vec, n);
//            mm_x_ = minmod(grid.sec_dx(ans_vec, n), grid.sec_dx(ans_vec, n-1));
//            dy_ = grid.dx_
//        }

//        if ( vx > 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_backward(ans_vec, n) + ( minmod( grid.sec_dx(ans_vec,n), grid.sec_dx(ans_vec,n-1) ) )*0.5*dx))  +
//                                                                            (vy * ( grid.dy_backward(ans_vec, n) + ( minmod( grid.sec_dy(ans_vec,n), grid.sec_dy(ans_vec,n-1) ) )*0.5*dy)) ); }
//        else if ( vx < 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_forward(ans_vec, n) + ( minmod( grid.sec_dx(ans_vec,n+1), grid.sec_dx(ans_vec,n) ) )*0.5*dx))  +
//                                                                                 (vy * ( grid.dy_forward(ans_vec, n) + ( minmod( grid.sec_dy(ans_vec,n+1), grid.sec_dy(ans_vec,n) ) )*0.5*dy)) ); }
//        else if ( vx < 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_forward(ans_vec, n) + ( minmod( grid.sec_dx(ans_vec,n+1), grid.sec_dx(ans_vec,n) ) )*0.5*dx))  +
//                                                                                 (vy * ( grid.dy_backward(ans_vec, n) + ( minmod( grid.sec_dy(ans_vec,n), grid.sec_dy(ans_vec,n-1) ) )*0.5*dy)) ); }
//        else if ( vx > 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_backward(ans_vec, n) + ( minmod( grid.sec_dx(ans_vec,n), grid.sec_dx(ans_vec,n-1) ) )*0.5*dx))  +
//                                                                                 (vy * ( grid.dy_forward(ans_vec, n) + ( minmod( grid.sec_dy(ans_vec,n+1), grid.sec_dy(ans_vec,n) ) )*0.5*dy)) ); }
//        else {ans_vec[n] = ans_vec[n];}

//        std::cout<<n <<std::endl;
        if ( vx >= 0.0 && vy >= 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_backward(ans_vec, n) + (( minmod( grid.sec_dx(ans_vec,n), grid.sec_dx(ans_vec,n-1) ) )*0.5*dx)))  +
                                                                            (vy * ( grid.dy_backward(ans_vec, n) + (( minmod( grid.sec_dy(ans_vec,n), grid.sec_dy(ans_vec,n-N) ) )*0.5*dy))) ); }
        else if ( vx < 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_forward(ans_vec, n) + (( minmod( grid.sec_dx(ans_vec,n+1), grid.sec_dx(ans_vec,n) ) )*0.5*dx)))  +
                                                                                 (vy * ( grid.dy_forward(ans_vec, n) + (( minmod( grid.sec_dy(ans_vec,n+N), grid.sec_dy(ans_vec,n) ) )*0.5*dy))) ); }
        else if ( vx < 0.0 && vy >= 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_forward(ans_vec, n) + (( minmod( grid.sec_dx(ans_vec,n+1), grid.sec_dx(ans_vec,n) ) )*0.5*dx)))  +
                                                                                 (vy * ( grid.dy_backward(ans_vec, n) + (( minmod( grid.sec_dy(ans_vec,n), grid.sec_dy(ans_vec,n-N) ) )*0.5*dy))) ); }
        else if ( vx >= 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] - dt*( (vx * ( grid.dx_backward(ans_vec, n) + (( minmod( grid.sec_dx(ans_vec,n), grid.sec_dx(ans_vec,n-1) ) )*0.5*dx)))  +
                                                                                 (vy * ( grid.dy_forward(ans_vec, n) + (( minmod( grid.sec_dy(ans_vec,n+N), grid.sec_dy(ans_vec,n) ) )*0.5*dy))) ); }
        else {ans_vec[n] = ans_vec[n];}

    }
}

void ENO_Advection::firstord(double dt)
{
    for (int n=0; n < (grid.get_N()*grid.get_M()) ; ++n)
    {
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);

        double dx = grid.get_dx();
        double dy = grid.get_dy();

        double vx = velocityx->operator ()(x,y);
        double vy = velocityy->operator ()(x,y);
        //double vx = (* velocityx)(x,y);
        //double vy = (* velocityy)(x,y);

        if ( vx > 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_backward(ans_vec, n) ))  +
                                                                            (vy * ( grid.dy_backward(ans_vec, n) )) ); }
        else if ( vx < 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_forward(ans_vec, n) ))  +
                                                                                 (vy * ( grid.dy_forward(ans_vec, n))) ); }
        else if ( vx < 0.0 && vy > 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_forward(ans_vec, n) ))  +
                                                                                 (vy * ( grid.dy_backward(ans_vec, n) )) ); }
        else if ( vx > 0.0 && vy < 0.0 ) { ans_vec[n] = ans_vec[n] + dt*(-1.0)*( (vx * ( grid.dx_backward(ans_vec, n) ))  +
                                                                                 (vy * ( grid.dy_forward(ans_vec, n) )) ); }
        else {ans_vec[n] = ans_vec[n];}


    }
}

