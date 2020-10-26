#include <iostream>
#include <math.h>
#include <grid2d.h>
#include <sl_method.h>
#include <cmath>
#include <omp.h>
#include <math_tools.h>
#include <algorithm>

int main()
{

    velocity_X x_velocity;
    velocity_Y y_velocity;

    int N = 128;
    Grid2D gridf(N,N,-1.0,1.0,-1.0,1.0);

    //double deltatau = 0.5*gridf.get_dx();
    double deltat = 0.5*gridf.get_dx();


    std::vector<double> insol;
    insol.resize(N*N);
    std::fill (insol.begin(),insol.end(),0.0);


    // INITIAL SOLUTION
    for (int n=0; n < (N*N) ; n++)
    {
        double x = gridf.x_from_n(n);
        double y = gridf.y_from_n(n);
        insol[n] = sqrt( (x - 0.25)*(x - 0.25) + (y*y) ) - 0.2;
    }

    //HEART
    for (int n=0; n < (N*N) ; n++)
    {
        double x = gridf.x_from_n(n);
        double y = gridf.y_from_n(n);
        if ( (x*x + 2.*( (3./5.)*pow(x*x, 1./3.) - y)*( (3./5.)*pow(x*x, 1./3.) - y) -1.) <= -0.4) {
            insol[n] = 0.;
        }
        else {insol[n] = 1. ;}
    }



//    // NUMBER FOUR: EXTRA CREDIT
//    for (int k=0; k < (150) ; k++) {
//        for (int n=0; n < (N*N) ; n++) {
//            insol[n] += k*0.001;
//        }
//        sprintf(name,"/Users/mpowell2/Documents/Home3/ghf22=%d.vtk",k);
//        gridf.initialize_VTK_file(name);
//        gridf.print_VTK_Format(insol, "value_at_nodes",name);
//    }

//    //PERTURBED SOLUTION
//    std::vector<double> insol_p;
//    insol_p.resize(N*N);
//    std::fill (insol_p.begin(),insol_p.end(),0.0);

//    for (int n=0; n < (N*N) ; ++n)
//    {
//        double x = gridf.x_from_n(n);
//        double y = gridf.y_from_n(n);
//        insol_p[n] = 5.*(sqrt( (x - 0.25)*(x - 0.25) + (y*y) ) - 0.2);
//    }



    //SL_method SL(gridf, insol_p, &x_velocity, &y_velocity);

    SL_method SL(gridf, insol, &x_velocity, &y_velocity);

    double tfinal = 2.0*M_PI;
    int num_iters = tfinal/deltat;
    int countmod = floor(num_iters / 50.);

    std::cout << num_iters << std::endl;

    std::vector <double> true_soln;
    true_soln.resize(N*N);
    std::fill (true_soln.begin(), true_soln.end(), 0.0);

    for(int n=0 ; n < (N*N) ; n++) {
        double x = gridf.x_from_n(n);
        double y = gridf.y_from_n(n);
        double t = tfinal;
        true_soln[n] = sqrt( ((x * cos(t) + y * sin(t)) - 0.25) * ((x * cos(t) + y * sin(t)) - 0.25)
                             + (y * cos(t) - x * sin(t)) * (y * cos(t) - x * sin(t)) ) - 0.2;
    }


//    NUMBER THREE
//    std::vector<double> sol(N*N);

//    for (int n=0; n < (num_iters+1) ; ++n)
//    {
//        sol = SL.get_sol();
//        SL.set_init(sol);
//        SL.one_step(deltat);
//        for (int n=0; n < (10) ; ++n ) {
//            SL.reinit(deltat/50.);
//            //std::vector<double> sol = SL.get_sol();
//        }
//        sol = SL.get_sol();
//        if (n  % countmod == 0 ) {
//            //char name[250];
//            sprintf(name,"/Users/mpowell2/Documents/FIN/finalthree2=%d.vtk",n);
//            gridf.initialize_VTK_file(name);
//            gridf.print_VTK_Format(sol, "value_at_nodes",name);
//        }
//    }


//    // NUMBER ONE
//    std::vector<double> sol(N*N);

//    for (int n=0; n < (num_iters) ; n++)
//    {
//        SL.one_step(deltat);
//        sol = SL.get_sol();
//        if (n  % countmod == 0 ) {
//            char name[250];
//            sprintf(name,"/Users/mpowell2/Documents/FIN/done2=%d.vtk",n);
//            gridf.initialize_VTK_file(name);
//            gridf.print_VTK_Format(sol, "value_at_nodes",name);
//        }
//    }


//    //CHECKING ERROR
//    sol = SL.get_sol();

//    std::vector<double> difference_vec;
//    difference_vec.resize(N*N);
//    std::fill (difference_vec.begin(),difference_vec.end(),0.0);

//    // TWO NORM
//    for (int p = 0; p<(N*N); p++)
//    {
//        difference_vec[p] = abs(sol[p] - true_soln[p]);
//    }

//    // MAX NORM
//    double sum;
//    for (int k=0; k< (N*N) ; k++)
//    {
//        sum = sum + (difference_vec[k]*difference_vec[k]);
//    }
//    double max_val = *max_element(difference_vec.begin(), difference_vec.end());

//    std::cout<< "L2: " << sqrt(sum) << std::endl;
//    std::cout<< "Max: " << max_val << std::endl;






//    //NUMBER TWO: REINITIALIZATION
//    std::vector<double> sol(N*N);
//    sol.resize(N*N);
//    for (int n=0; n < (250) ; ++n)
//    {
//        SL.reinit(deltatau);
//        sol = SL.get_sol();
//        sprintf(name,"/Users/mpowell2/Documents/FIN/reinit1=%d.vtk",n);
//        gridf.initialize_VTK_file(name);
//        gridf.print_VTK_Format(sol, "value_at_nodes",name);
//    }


    std::cout << "Hello World!" << std::endl;
    return 0.;

}
