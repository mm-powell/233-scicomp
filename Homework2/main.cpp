#include <iostream>
#include <cf_2.h>
#include <math.h>
#include <grid2d.h>
#include <eno_advection.h>
#include <cmath>
#include <omp.h>

using namespace std;

//class velocity_X :CF_2
//{
//public:
//    double operator()(double x, double y) const{
//        return cos(x)*sin(y);
//    }
//};


//class velocity_Y :CF_2
//{
//public:
//    double operator()(double x, double y) const{
//        return 1.;
//    }
//};

//double vec_norm(std::vector<double> diff_vec);

//std::vector<double> diff(std::vector<double> approx_vec, std::vector<double> true_vec);

//double vec_norm(std::vector<double> diff_vec) {
//    // two norm
//    // square each element, then take square root of entire thing
//    double sum;
//    for (int p=0; p<=int(diff_vec.size()); ++p)
//    {
//        sum = sum + (diff_vec[p]*diff_vec[p]);
//    }
//    return sqrt(sum);
//}


//std::vector<double> diff(std::vector<double> approx_vec, std::vector<double> true_vec) {
//    // difference between error vector and true solution
//    std::vector<double> difference_vec;
//    difference_vec.resize(approx_vec.size());
//    std::fill (difference_vec.begin(),difference_vec.end(),0.0);

//    for (int p = 0; p<=(approx_vec.size()); ++p)
//    {
//        difference_vec[p] = std::abs(approx_vec[p] - true_vec[p]);
//    }

//    return difference_vec;
//}


int main()
{

    // creating velocities (pointers)
    velocity_X x_velocity;
    velocity_Y y_velocity;

    // creating a grid
    int N = 64;
    Grid2D gridf(N,N,-1.0,1.0,-1.0,1.0);

    // dt
    //double deltat = (1./64.)*gridf.get_dx();
    double deltat = 0.25*gridf.get_dx();
    //double deltat = (1./32.) * gridf.get_dx();
    //double deltat = gridf.get_dx()*gridf.get_dx();

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

    //int countmod = floor(num_iters / 300.);

    cout << tfinal << endl;

    cout << num_iters << endl;

    for (int n=0; n < (num_iters + 2) ; ++n)
    {
        EA.Central_Advection(deltat);
        std::vector<double> solution = EA.get_sol();

        //if (n % countmod == 0) {
            char name[250];
            sprintf(name,"/Users/mpowell2/Documents/CA/N=%d.vtk",n);
            gridf.initialize_VTK_file(name);
            gridf.print_VTK_Format(solution, "value_at_nodes",name);
        //}
        cout << n << endl;
    }





//    std::vector <double> true_soln;
//    true_soln.resize(N*N);
//    std::fill (true_soln.begin(), true_soln.end(), 0.0);

//    for(int n=0 ; n < (N*N) ; ++n) {
//        double x = gridf.x_from_n(n);
//        double y = gridf.y_from_n(n);
//        double t = tfinal;
//        //double t = n*deltat;
//        if ( sqrt( ( (x*cos(t) + y*sin(t) - 0.5)*(x*cos(t) + y*sin(t) - 0.5) ) + ( (y*cos(t) - x*sin(t))*(y*cos(t) - x*sin(t)) ) ) <= 0.2 ) {true_soln[n] = 1.0;}
//        else {true_soln[n] = 0.0;}
////        char name[250];
////        sprintf(name,"/Users/mpowell2/Documents/true/true_%d.vtk", n);
////        gridf.initialize_VTK_file(name);
////        gridf.print_VTK_Format(true_soln, "value_at_nodes",name);
//    }


    //std::vector<double> ans_vec_2 = EA.get_sol();

//    std::vector<double> difference_vec;
//    difference_vec.resize(N*N);
//    std::fill (difference_vec.begin(),difference_vec.end(),0.0);

//    for (int p = 0; p<(N*N); ++p)
//    {
//        //double placeh = ans_vec_2[p] - true_soln[p];
//        difference_vec[p] = std::abs(ans_vec_2[p] - true_soln[p]);
//    }

//    double sum;
//    for (int p=0; p< (N*N) ; ++p)
//    {
//        sum = sum + (difference_vec[p]*difference_vec[p]);
//    }

//    cout<< sqrt(sum) <<endl;



//    for (int n=0; n < (64*64) ; ++n)
//    {
//        cout<<ans[n]<<endl;
//    }


    //TIME ORDER
//    std::vector<double> deltat_vec{( (1./4.) * gridf.get_dx() ), ( (1./8.) * gridf.get_dx() ), ( (1./16.) * gridf.get_dx() ), ( (1./32.) * gridf.get_dx() )};

//    std::vector<double> error_t_vec;
//    error_t_vec.resize(deltat_vec.size());
//    for(int k=0 ; k < deltat_vec.size() ; ++k) {
//        ENO_Advection EA(gridf, ans, &x_velocity, &y_velocity);
//        for (int n=0; n < (num_iters + 2) ; ++n)
//        {
//            EA.one_Step(deltat_vec[k]);
//        }
//        std::vector<double> solution = EA.get_sol();
//        error_t_vec[k] = vec_norm( diff(solution, true_soln) );
//        cout<<error_t_vec[k]<<endl;
//    }







    //SPACE ORDER



    cout << "Hello World!" << endl;
    return 0;
}




