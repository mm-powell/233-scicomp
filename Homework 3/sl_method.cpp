#include "sl_method.h"
#include <grid2d.h>
#include <math_tools.h>

SL_method::SL_method()
{

}


SL_method::SL_method(Grid2D grid_, std::vector<double> init_sol_, velocity_X *velocity_x_, velocity_Y *velocity_y_){
    grid = grid_;
    solution.resize(grid.get_N()*grid.get_M());
    solution = init_sol_;
    velocity_x = velocity_x_;
    velocity_y = velocity_y_;
    init_sol.resize(grid.get_N()*grid.get_M());
    init_sol = init_sol_;
}


void SL_method::set_init(std::vector<double> sol)
{
    init_sol = sol;
}


std::vector<double> SL_method::trajectory_interpolation(int n, double dt){

    double xmin = grid.get_xmin();
    double xmax = grid.get_xmax();

    double ymin = grid.get_ymin();
    double ymax = grid.get_ymax();

    double x_n = grid.x_from_n(n);
    double y_n = grid.y_from_n(n);

    double x_star = ( x_n ) - (dt/2.)*( (* velocity_x)(x_n,y_n) );
    double y_star = ( y_n ) - (dt/2.)*( (* velocity_y)(x_n,y_n) );

    if (x_star <= xmin) { x_star = xmin; }
    if (x_star >= xmax ) { x_star = xmax; }

    if (y_star <= ymin) { y_star = ymin; }
    if (y_star >= ymax ) { y_star = ymax; }

    double x_d = ( x_n ) - (dt)*( (* velocity_x)( x_star , y_star) );
    double y_d = ( y_n ) - (dt)*( (* velocity_y)( x_star , y_star) );

    if (x_d <= xmin) { x_d = xmin; }
    if (x_d >= xmax ) { x_d = xmax; }

    if (y_d <= ymin) { y_d = ymin; }
    if (y_d >= ymax ) { y_d = ymax; }

    std::vector<double> coord;
    coord.resize(2);

    coord[0] = x_d;
    coord[1] = y_d;

    return coord;
}

void SL_method::one_step(double dt){

//#pragma omp parallel for
    for (int n=0; n < (grid.get_N()*grid.get_M()) ; n++) {
        std::vector<double> depart = trajectory_interpolation( n, dt );
        double x_d = depart[0];
        double y_d = depart[1];

        solution[n] = ENO_interpolation(grid, solution, x_d, y_d);
        //init_sol[n] = solution[n];

    }

}

std::vector<double> SL_method::get_sol()
{
    return solution;
}

void SL_method::reinit(double dt) {
    std::vector<double> func = solution;
//#pragma omp parallel for
    for (int n=0; n < (grid.get_N()*grid.get_M()) ; n++)
    {
        double phi_x = 0.;
        double phi_y = 0.;
        //double phi_x;
        //double phi_y;

        if ( grid.signum(init_sol[n])*grid.dx_backward(func, n) <= 0.0 && grid.signum(init_sol[n])*grid.dx_forward(func, n) <= 0.0 )
        { phi_x = grid.dx_forward(func, n); }
        else if ( grid.signum(init_sol[n])*grid.dx_backward(func, n) >= 0.0 && grid.signum(init_sol[n])*grid.dx_forward(func, n) >= 0.0 )
        { phi_x = grid.dx_backward(func, n); }
        else if ( grid.signum(init_sol[n])*grid.dx_backward(func, n) <= 0.0 && grid.signum(init_sol[n])*grid.dx_forward(func, n) >= 0.0 )
        { phi_x = 0.; }
        //else
        else if ( grid.signum(init_sol[n])*grid.dx_backward(func, n) >= 0.0 && grid.signum(init_sol[n])*grid.dx_forward(func, n) <= 0.0 )
        {
            if ( abs(grid.dx_backward(func, n)) >= abs(grid.dx_forward(func, n)) ) { phi_x = grid.dx_backward(func, n); }
            //else { phi_x = grid.dx_forward(func, n); }
            else if ( abs(grid.dx_backward(func, n)) <= abs(grid.dx_forward(func, n)) ) { phi_x = grid.dx_forward(func, n); }
        }

        if ( grid.signum(init_sol[n])*grid.dy_backward(func, n) <= 0.0 && grid.signum(init_sol[n])*grid.dy_forward(func, n) <= 0.0 )
        { phi_y = grid.dy_forward(func, n); }
        else if ( grid.signum(init_sol[n])*grid.dy_backward(func, n) >= 0.0 && grid.signum(init_sol[n])*grid.dy_forward(func, n) >= 0.0 )
        { phi_y = grid.dy_backward(func, n); }
        else if ( grid.signum(init_sol[n])*grid.dy_backward(func, n) <= 0.0 && grid.signum(init_sol[n])*grid.dy_forward(func, n) >= 0.0 )
        { phi_y = 0.; }
        //else
        else if ( grid.signum(init_sol[n])*grid.dy_backward(func, n) >= 0.0 && grid.signum(init_sol[n])*grid.dy_forward(func, n) <= 0.0 )
        {
            if ( abs(grid.dy_backward(func, n)) >= abs(grid.dy_forward(func, n)) ) { phi_y = grid.dy_backward(func, n); }
            //else { phi_y = grid.dy_forward(func, n); }
            else if ( abs(grid.dy_backward(func, n)) <= abs(grid.dy_forward(func, n)) ) { phi_y = grid.dy_forward(func, n); }
        }


        solution[n] = func[n] - dt*grid.signum(init_sol[n])*( sqrt( phi_x*phi_x  + phi_y*phi_y ) - 1. );

    }
}


