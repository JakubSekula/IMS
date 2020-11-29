#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double b = 1.2;
const double g = 0.04;

double N = 1000.0;

double eta = 0.5;
double sigma = 1/5.1;
double alpha = 0.5;
double beta = 0.5;
double greeky = 0.025;

double iotaA = 1/7.0;
double iotaI = 1/7.0;
double iotaH = 1/14.0;

double delta = 0.015;
double epsilon = 1.0;

//typedef boost::array< double , 6 > state_type;
typedef boost::array< double , 14 > state_type;

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

    N = x[ 0 ] + x[ 1 ] + x[ 2 ] + x[ 3 ] + x[ 5 ];

}

void sir( const state_type &x , state_type &dxdt , double t )
{
    t = t;
    dxdt[0] = ( -b *  x[0] * x[1] ) / N;
    dxdt[1] = ( b *  ( x[0] * x[1] ) / N )  - g * x[1];
    dxdt[2] = g * x[1];
}

void write_sir(const state_type &x , const double t )
{
    cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[ 4 ] << ' ' << x[ 5 ] << ' ' << x[ 6 ] << endl;
}

void write_sir_masked(const state_type &x , const double t )
{
    cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[ 4 ] << ' ' << x[ 5 ] << ' ' << x[ 6 ] << ' ' << x[7] << ' ' << x[8] << ' ' << x[9] << ' ' << x[10] << ' ' << x[ 11 ] << ' ' << x[ 12 ] << ' ' << x[ 13 ] << endl;
}

void masked_real( const state_type &y, state_type &dxdt, double t ){
    t = t;
    dxdt[0] = -b * ( y[2] + eta * y[3] ) * ( y[0] / N ) - b * ( ( 1 - epsilon ) * y[10] + ( 1 - epsilon ) * eta * y[8] ) * ( y[0] / N );
    dxdt[1] = b * ( y[2] + eta * y[3] ) * ( y[0] / N ) - b * ( ( 1 - epsilon ) * y[10] + ( 1 - epsilon ) * eta * y[8] ) * ( y[0] / N ) - ( sigma * y[1] );
    dxdt[2] = alpha * sigma * y[1] - greeky * y[2] - iotaI * y[2];
    dxdt[3] = ( 1 - alpha ) * sigma * y[1] - iotaA * y[3];
    dxdt[4] = greeky * y[2] - delta * y[4] - iotaH * y[4];
    dxdt[5] = iotaI * y[2] + iotaA * y[3] + iotaH * y[4];
    dxdt[6] = delta * y[4];

    dxdt[7] = -b * ( 1 - epsilon ) * ( y[2] + eta * y[3] ) * ( y[7] / N ) - b * ( 1 - epsilon ) * ( ( 1 - epsilon ) * y[9] + ( 1 - epsilon ) * eta * y[10] ) * ( y[7] / N );
    dxdt[8] = b * ( 1 - epsilon ) * ( y[2] + eta * y[3] ) * ( y[7] / N ) + b * ( 1 - epsilon ) * ( ( 1 - epsilon ) * y[9] + ( 1 - epsilon ) * eta * y[10] ) * ( y[7] / N ) - sigma * y[8];
    dxdt[9] = alpha * sigma * y[8] - greeky * y[9] - iotaI * y[9];
    dxdt[9] = greeky * y[9] - delta * y[11] - iotaH * y[11];
    dxdt[10] = iotaI * y[9] + iotaA * y[10] + iotaH * y[11];
    dxdt[11] = delta * y[11];

    N = y[ 0 ] + y[ 1 ] + y[ 2 ] + y[ 3 ] + y[ 5 ] + y[ 7 ] + y[ 8 ] + y[ 9 ] + y[ 10 ] + y[ 12 ];

}

int main(int argc, char **argv)
{
    argc = argc;
    argv = argv;
    //                 S     E    I    A    H    R    D
    //state_type x = { N - 6, 2.0, 3.0, 1.0, 0.0, 1.1, 0.0 }; // initial conditions
    //integrate( masked , x , 0.0 , 200.0 , 0.1 , write_sir );

    //               Su[0]  Eu[1]   Iu[2]   Au[3]   Hu[4]   Ru[5]   Du[6]   Sm[7]   Em[8]   Im[9]   Am[10]   Hm[11]  Rm[12] Dm[13]
    state_type y = { N - 6,   2.0,    3.0,    1.0,    0.0,    1.1,    0.0,    N - 6,   2.0,    3.0,    1.0,    0.0,    1.1,    0.0  };
    integrate( masked_real, y, 0.0, 200.0, 0.1, write_sir_masked );

}