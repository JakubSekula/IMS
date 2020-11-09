#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

double b;
const double g = 0.05;

typedef boost::array< double , 3 > state_type;



void getArgs( int argc, char** argv ){

    int arg;
    string temp = "";
    while( ( arg = getopt( argc, argv, "r:" ) ) != -1 ){
        switch( arg ){
            case 'r':
                b = stod( optarg );
                break;
            default:
                fprintf( stderr, "Uknown argument\n" );
                exit( 10 );
        }
    }

    if( !b ){
        fprintf( stderr, "Reproduction number must be entered\n" );
        exit( 10 );
    }

}

void sir( const state_type &x , state_type &dxdt , double t )
{
    t = t; 
    dxdt[0] = -b * x[0] * x[1];
    dxdt[1] = b * x[0] * x[1] - g * x[1];
    dxdt[2] = g * x[1];
}

void write_sir(const state_type &x , const double t )
{
    cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
}

int main( int argc, char **argv )
{
    getArgs( argc, argv );

    state_type x = { 0.99 , 0.01 , 0.0 }; // initial conditions
    integrate( sir , x , 0.0 , 200.0 , 0.1 , write_sir );
}