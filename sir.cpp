#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double b = 0.4;
const double g = 0.04;

double N = 1000.0;

double eta = 0.5;
double sigma = 1/5.1;
double alpha = 0.5;
double beta = 0.5;
double greeky = 0.025;

double iotaA = 1/7;
double iotaI = 1/7;
double iotaH = 1/14;

double delta = 0.015;

typedef boost::array< double , 7 > state_type;

void masked( const state_type &x , state_type &dxdt , double t )
{
    t = t; 
    dxdt[0] = -b * ( ( x[ 2 ] + eta * x[ 3 ] ) * ( x[ 0 ] / N ) ) ;
    dxdt[1] = b * ( ( x[ 2 ] + eta * x[ 3 ] ) * ( x[ 0 ] / N ) ) - sigma * x[ 1 ];
    dxdt[2] = alpha * sigma * x[ 1 ] - greeky * x[ 2 ] - iotaI * x[ 2 ];
    dxdt[3] = ( 1 - alpha ) * sigma * x[ 1 ] - iotaA * x[ 3 ];
    dxdt[4] = greeky * x[ 2 ] - delta * x[ 4 ] - iotaH * x[ 4 ];
    dxdt[5] = iotaI * x[ 2 ] + iotaA * x[ 3 ] + iotaH * x[ 4 ];
    dxdt[6] = delta * x[ 4 ];
}

void sir( const state_type &x , state_type &dxdt , double t )
{
    t = t;
    dxdt[0] = ( -b *  x[0] * x[1] ) / N;
    dxdt[1] = ( b *  ( x[0] * x[1] ) / N )  - g * x[1];
    dxdt[2] = g * x[1];

    N = x[ 0 ] + x[ 1 ] + x[ 2 ] + x[ 3 ] + x[ 5 ];

}

void write_sir(const state_type &x , const double t )
{
    cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[ 4 ] << ' ' << x[ 5 ] << ' ' << x[ 6 ] << endl;
}

int main(int argc, char **argv)
{
    argc = argc;
    argv = argv;
    //                 S     E    I    A    H    R    D
    state_type x = { N - 6, 2.0, 3.0, 1.0, 0.0, 1.0, 0.0 }; // initial conditions
    integrate( masked , x , 0.0 , 200.0 , 0.1 , write_sir );
}