#include <iostream>
#include <vector>

using namespace std;

double Legendre ( double x, int n )
{
    if ( n == 1 ) { return x; }
    else if ( n == 2 ) { return ( 1.0 / 2.0 ) * (( 3.0 * x * x) - 1.0); }
    else if ( n == 3 ) { return ( 1.0 / 2.0 ) * ((5.0 * x * x * x) - (3.0 * x));}
    else if ( n == 4 ) { return ( 1.0 / 8.0 ) * ((35.0 * x * x * x * x) - (30.0 * x * x) + 3);}
    else if ( n == 5 ) { return ( 1.0 / 8.0 ) * ((63.0 * x * x * x * x * x) - (70.0 * x * x * x) + ( 15.0 * x ));}
    else  { cout<<"Error!"<<endl; return 1; }
}

std::vector<double> sampledLegendre ( double a, double b, int N, int n )
{
    // create a vector of size N to store answers
    std::vector<double> ans_vec;
    ans_vec.resize(N);

    //define step size based on a, b, and N
    double step = (b-a)/(N-1);

    // create a vector of size N to fill with uniformly spaced x values
    std::vector<double> X;
    X.resize(N);

    // populate X vector with uniformly space x values
    for (int i=0; i<N; ++i)
    {
        X.at(i)= (a + (step * i));
    }

    for (int j=0; j<N; ++j)
    {
        cout<<X[j]<<endl;
    }

    // populating answer vector with Legendre values by calling previous function
    for (int k = 0 ; k < N ; ++k)
    {
        ans_vec.at(k) = Legendre(X[k], n);
    }

    //for (int p=0; p<N; ++p)
    //{
    //    cout<<ans_vec[p]<<endl;
    //}

    return ans_vec;
}


int main()
{

    //double ans = Legendre(0,6);
    //cout<<ans<<endl;

    std::vector<double> ans_vec_;
    ans_vec_  = sampledLegendre(0,2,5,8);

    for (int p=0; p<int(ans_vec_.size()); ++p)
    {
        cout<<ans_vec_[p]<<endl;
    }

    return 0;
}
