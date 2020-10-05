#include <iostream>
#include <cf_2.h>
#include <math.h>
#include <grid2d.h>
#include <eno_advection.h>
#include <cmath>
#include <omp.h>

using namespace std;


int main()
{

    // creating velocities (pointers)
    velocity_X x_velocity;
    velocity_Y y_velocity;

    // creating a grid
    int N = 64;
    Grid2D gridf(N,N,-1.0,1.0,-1.0,1.0);

    // dt
    double deltat = 0.25*gridf.get_dx();

    // creating empty answer vector
    std::vector<double> ans;
    ans.resize(N*N);
    std::fill (ans.begin(),ans.end(),0.0);

    // initial solution
    for (int n=0; n < (N*N) ; ++n)
    {
        double x = gridf.x_from_n(n);
        double y = gridf.y_from_n(n);
        if ( ( sqrt( ((x - (0.5)) * (x - (0.5))) + (y*y) ) - 0.2 ) <= 0.0 ) { ans[n] = 1.0; }
        else { ans[n] = 0.0; }
        //ans[n] = sqrt( (x - 0.25)*(x - 0.25) + (y*y) ) - 0.2;
    }


    ENO_Advection EA(gridf, ans, &x_velocity, &y_velocity);

    double tfinal = 2.0*M_PI;
    int num_iters = tfinal/deltat;

    cout << tfinal << endl;
    cout << num_iters << endl;

    // primary loop
    for (int n=0; n < (num_iters + 2) ; ++n)
    {
        EA.Central_Advection(deltat);
        std::vector<double> solution = EA.get_sol();

        char name[250];
        sprintf(name,"/Users/mpowell2/Documents/CA/N=%d.vtk",n);
        gridf.initialize_VTK_file(name);
        gridf.print_VTK_Format(solution, "value_at_nodes",name);
        cout << n << endl;
    }


    std::vector <double> true_soln;
    true_soln.resize(N*N);
    std::fill (true_soln.begin(), true_soln.end(), 0.0);

    for(int n=0 ; n < (N*N) ; ++n) {
        double x = gridf.x_from_n(n);
        double y = gridf.y_from_n(n);
        double t = tfinal;
        if ( sqrt( ( (x*cos(t) + y*sin(t) - 0.5)*(x*cos(t) + y*sin(t) - 0.5) ) + ( (y*cos(t) - x*sin(t))*(y*cos(t) - x*sin(t)) ) ) <= 0.2 ) {true_soln[n] = 1.0;}
        else {true_soln[n] = 0.0;}
//        char name[250];
//        sprintf(name,"/Users/mpowell2/Documents/true/true_%d.vtk", n);
//        gridf.initialize_VTK_file(name);
//        gridf.print_VTK_Format(true_soln, "value_at_nodes",name);
    }


    std::vector<double> difference_vec;
    difference_vec.resize(N*N);
    std::fill (difference_vec.begin(),difference_vec.end(),0.0);

    for (int p = 0; p<(N*N); ++p)
    {
        difference_vec[p] = std::abs(ans[p] - true_soln[p]);
    }

    double sum;
    for (int p=0; p< (N*N) ; ++p)
    {
        sum = sum + (difference_vec[p]*difference_vec[p]);
    }


    cout << "Hello World!" << endl;
    return 0;
}




