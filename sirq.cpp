#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double b = 0.13; // No of people falling ill
const double g = 0.05; // No of peple recovering
const double ft = 0.2; // Fraction of infected people tested

typedef boost::array<double, 4> state_type;
typedef boost::array<double, 6> state_type2;

void sir(const state_type &x , state_type &dxdt , double t)
{
    dxdt[0] = -b * ((x[0] * x[1]) / x[3]);
    dxdt[1] = b * ((x[0] * x[1]) / x[3]) - g * x[1];
    dxdt[2] = g * x[1];
}

void sirq(const state_type2 &x, state_type &dxdt, double t)
{
    dxdt[0] = -b * ((x[0] * x[4] * x[3]) / x[5]);
    dxdt[1] = ft * (b * ((x[0] * x[4] * x[3]) / x[5])) - g * x[3];
    dxdt[2] = (1 - ft) * (b * ((x[0] * x[4] * x[3]) / x[5])) - g * x[3] * x[4];
    dxdt[3] = g * x[3] + g * x[3] * x[4];
}

void write_sir(const state_type &x , const double t )
{
    cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
}

int main(int argc, char **argv)
{
    //state_type x = { 4999999 , 1 , 0 , 5000000 }; // initial conditions
    //integrate( sir , x , 0.0 , 500.0 , 0.1 , write_sir );

    //                   S      I   R   Q      U         N
    state_type2 y = { 4999999 , 1 , 0 , 0 , 5000000 , 5000000 };
    integrate(sirq, y, 0.0, 500.0, 0.1, write_sir);

}
