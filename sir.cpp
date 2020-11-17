#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double b = 0.4236;
const double g = 0.04;

double N = 14368332;

double alpha = 0.5;
double ny = 0.5;
double sigma = 1/5.2;
double theta = 1.8999 * pow( 10, -12 );
double psi = 0.0135;

double delta = 0.0;
double epsilon = 0.0;

double iotaA = 0.13978;
double iotaO = 0.13978;
double iotaI = 1/15;

double d0 = 0.1113;

typedef boost::array< double , 6 > state_type;
//typedef boost::array< double , 14 > state_type;

void sir( const state_type &x , state_type &dxdt , double t )
{
    t = t;
    dxdt[0] = ( -b *  x[0] * x[1] ) / N;
    dxdt[1] = ( b *  ( x[0] * x[1] ) / N )  - g * x[1];
    dxdt[2] = g * x[1];
}

void write_sir(const state_type &x , const double t )
{
    cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[ 4 ] << ' ' << x[ 5 ] << endl;
}

void masked( const state_type &x , state_type &dxdt , double t )
{
    t = t; 
    dxdt[0] = - ( ( b * ( 1 - delta ) * ( 1 - epsilon ) * ( alpha * x[2] + x[3] ) ) / ( N - x[4] ) ) * x[0];
    dxdt[1] = ( ( b * ( 1 - delta ) * ( 1 - epsilon ) * ( alpha * x[2] + x[3] ) ) / ( N - x[4] ) ) * x[0] - sigma * x[1];
    dxdt[2] = ny * sigma * x[1] - ( theta + iotaA ) * x[2];
    dxdt[3] = ( 1 - ny ) * sigma * x[1] - ( psi + iotaO + d0 ) * x[3];
    dxdt[4] = theta * x[2] + psi * x[3] - ( iotaI + d0 ) * x[4];
    dxdt[5] = iotaI * x[4] + iotaA * x[2] + iotaO * x[3];

    N = x[0] + x[1] + x[2] + x[3] + x[4] + x[5];

}

int main(int argc, char **argv)
{
    argc = argc;
    argv = argv;
    //                S[0]    E[1]    A[2]    I[3]       Id[4]   R[5]
    state_type x = {   N - 842, 441.0,   188.0, 212.0,     1.0,   0.0 }; // initial conditions
    integrate( masked , x , 0.0 , 400.0 , 0.1 , write_sir );

}