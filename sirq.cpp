#include <iostream>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

double b = 0.7; // No of people falling ill
double v = 0.362; // No of people recovering
double ft = 0.32; // Fraction of infected people tested
int population = 34000000; // Total population

// Defining a shorthand for the type of the mathematical state
typedef std::vector< double > state_type;

// SIR = suspectible, infected, recovered
void sir( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = (-b * x[0] * x[1]) / population;
    dxdt[1] = ((b * x[0] * x[1]) / population) - v * x[1];
    dxdt[2] = v * x[1];
}

// SIRQ = suspectible, infected, recovered, quarantined
void sirq( const state_type &x , state_type &dxdt , double t )
{    
    dxdt[0] = (-(b * x[0] * x[2])) / population;
    dxdt[1] = (ft * ((b * x[0] * x[2]) / population)) - (v * x[1]);
    dxdt[2] = ((1.0 - ft) * ((b * x[0] * x[2]) / population)) - (v * x[2]);
    dxdt[3] = (v * x[1]) + (v * x[2]);
}

// Observer, prints time and state when called (during integration)
void my_observer( const state_type &x, const double t )
{    
    //cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << endl;
    cout << t << ' ' << 0 << ' ' << x[1] + x[2] << ' ' << 0 << endl;
}


int main()
{
    // Initial state
    state_type x(4);
    x[0] = population - 2; // S
    x[1] = 0.0; // Q
    x[2] = 2.0; // UQ
    x[3] = 0.0; // R

    // Integration parameters
    double t0 = 0.0;
    double t1 = 18.0;
    double dt = 0.1;

    // Run integrator
    integrate(sirq, x, t0, t1, dt, my_observer);

    // Integration parameters
    t0 = 18.0;
    t1 = 55.0;
    dt = 0.1;

    b = 0.46; // No of people falling ill
    ft = 0.5; // Fraction of infected people tested

    // Run integrator
    integrate(sirq, x, t0, t1, dt, my_observer);

    // Integration parameters
    t0 = 55.0;
    t1 = 100.0;
    dt = 0.1;

    b = 0.7; // No of people falling ill
    ft = 0.32; // Fraction of infected people tested

    // Run integrator
    integrate(sirq, x, t0, t1, dt, my_observer);
}
