#ifndef SL_METHOD_H
#define SL_METHOD_H

#include <grid2d.h>
#include <vector>
#include <cf_2.h>


class SL_method
{
private:
    Grid2D grid;
    std::vector<double> solution;
    velocity_X * velocity_x;
    velocity_Y * velocity_y;
    std::vector<double> init_sol;

public:
    SL_method();
    SL_method(Grid2D grid_, std::vector<double> init_sol_, velocity_X *velocity_x_, velocity_Y *velocity_y_);
    std::vector<double> trajectory_interpolation(int n, double dt);
    void one_step(double dt);
    std::vector<double> get_sol();
    void reinit(double dt);
    void set_init(std::vector<double> sol);
};

#endif // SL_METHOD_H
