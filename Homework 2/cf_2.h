#ifndef CF_2_H
#define CF_2_H

#include <math.h>


class CF_2 //Continuous Function of 2 Variables
        //has contructor and operator that takes in two doubles and returns doubles
{
public:
    //CF_2();
    virtual double operator()(double x, double y) const = 0;
    // creating template class
    // Why use virtual? -- because we're going to overload later
};

class velocity_X :CF_2
{
public:
    double operator()(double x, double y) const{
        return -1.*y;
    }
};


class velocity_Y :CF_2
{
public:
    double operator()(double x, double y) const{
        return x;
    }
};

#endif // CF_2_H
