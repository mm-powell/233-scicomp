#ifndef ENO_ADVECTION_H
#define ENO_ADVECTION_H

#include <vector>
#include <grid2d.h>
#include <cf_2.h>

class ENO_Advection
{
private:
    Grid2D grid;
    std::vector<double> ans_vec;
    velocity_X * velocityx;
    velocity_Y * velocityy;
public:
    ENO_Advection();
    ENO_Advection(Grid2D grid_, std::vector<double> ans_vec_, velocity_X * velocityx_, velocity_Y * velocityy_);
    double minmod(double vel1, double vel2);
    void one_Step(double dt);
    void Central_Advection(double dt);
    void firstord(double dt);
    std::vector<double> get_sol();
    //ENO_Advection(double xmin, double xmax, double ymin, double ymax, double N);
    //std::vector<double> ans_vec;
    //double * x_direc;
    //double * y_direc;
};

#endif // ENO_ADVECTION_H
