#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// Transmission rate due to contacts with UNDETECTED asymptomatic infected
double alfa = 0.57;
// Transmission rate due to contacts with DETECTED asymptomatic infected
double beta1 = 0.0114;
// Transmission rate due to contacts with UNDETECTED symptomatic infected
double gamma1 = 0.456;
// Transmission rate due to contacts with DETECTED symptomatic infected
double delta = 0.0114;

// Detection rate for ASYMPTOMATIC
double epsilon = 0.171;
// Detection rate for SYMPTOMATIC
double theta = 0.3705;

// Worsening rate: UNDETECTED asymptomatic infected becomes symptomatic
double zeta = 0.1254;
// Worsening rate: DETECTED asymptomatic infected becomes symptomatic
double eta = 0.1254;

// Worsening rate: UNDETECTED symptomatic infected develop life-threatening symptoms
double mu = 0.0171;
// Worsening rate: DETECTED symptomatic infected develop life-threatening symptoms
double nu = 0.0274;

// Mortality rate for infected with life-threatening symptoms
double tau = 0.01;

// Recovery rate for undetected asymptomatic infected
double lambda = 0.0342;
// Recovery rate for detected asymptomatic infected
double rho = 0.0342;
// Recovery rate for undetected symptomatic infected
double kappa = 0.0171;
// Recovery rate for detected symptomatic infected
double xi = 0.0171;
// Recovery rate for life-threatened symptomatic infected
double sigma = 0.0171;

double r1 = epsilon + zeta + lambda;

double r2 = eta + rho + kappa;

double r3 = theta + mu + kappa;

double r4 = nu + xi;

double r5 = sigma + tau;


typedef boost::array<double, 8> state_type;

void sidarthe(const state_type &x , state_type &dxdt , double t)
{
    dxdt[0] = -x[0] * (alfa * x[1] + beta1 * x[2] + gamma1 * x[3] + delta * x[4]);
    dxdt[1] = x[0] * (alfa * x[1] + beta1 * x[2] + gamma1 * x[3] + delta * x[4]) - (epsilon + zeta + lambda) * x[1];
    dxdt[2] = epsilon * x[1] - (eta + rho) * x[2];
    dxdt[3] = zeta * x[1] - (theta + mu + kappa) * x[3];
    dxdt[4] = eta * x[2] + theta * x[3] - (nu + xi) * x[4];
    dxdt[5] = mu * x[3] + nu * x[4] - (sigma + tau) * x[5];
    dxdt[6] = lambda * x[1] + rho * x[2] + kappa * x[3] + xi * x[4] + sigma * x[5];
    dxdt[7] = tau * x[5];
}

void write_sidarthe(const state_type &x , const double t)
{
    cout << t << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[ 4 ] << ' ' << x[ 5 ] << endl;
    double R0 = ( alfa / r1 ) + ( ( beta1 * epsilon ) / ( r1 * r2 ) ) + ( ( gamma1 * zeta ) / ( r1 * r3 ) ) + ( ( delta * eta * epsilon ) / ( r1 * r2 * r4 ) ) + ( ( delta * zeta * theta ) / ( r1 * r3 * r4 ) );
    //cout << R0 << "\n";
}

int main(int argc, char **argv)
{
    int population = 60000000;
    double S0 = 1.0 - 200.0/population - 20.0/population - 1.0/population - 2.0/population - 0.0 - 0.0 - 0.0;
    
    // Initial state
    //             S[0]     I[1]              D[2]            A[3]           R[4]        T[5] H[6] E[7]
    state_type x = {S0, 200.0/population, 20.0/population, 1.0/population, 2.0/population, 0.0, 0.0, 0.0};
    

    
    alfa = 0.36;
    beta1 = 0.005;
    gamma1 = 0.2;
    delta = 0.005;
        
    mu = 0.008;
    nu = 0.015;
        
    zeta = 0.034;
    eta = 0.034;
    
    lambda = 0.08;
    rho = 0.0171;
    kappa = 0.0171;
    xi = 0.0171;
    sigma = 0.0171;
    // Days 0-4
    integrate(sidarthe, x, 0.0, 22.0, 0.01, write_sidarthe);

    /* // Days 4-12
    alfa = 0.4218;
    gamma1 = 0.285;
    beta1 = 0.0057;
    delta = 0.0057;
    integrate(sidarthe, x, 4.0, 12.0, 0.01, write_sidarthe);

    // Days 12-22
    epsilon = 0.1425;
    integrate(sidarthe, x, 12.0, 22.0, 0.01, write_sidarthe);

    // Days 22-28
    alfa = 0.36;
    beta1 = 0.005;
    gamma1 = 0.2;
    delta = 0.005;
        
    mu = 0.008;
    nu = 0.015;
        
    zeta = 0.034;
    eta = 0.034;
    
    lambda = 0.08;
    rho = 0.0171;
    kappa = 0.0171;
    xi = 0.0171;
    sigma = 0.0171;
    integrate(sidarthe, x, 22.0, 28.0, 0.01, write_sidarthe);

    // Days 28-38
    alfa = 0.21;
    gamma1 = 0.11;
    integrate(sidarthe, x, 28.0, 38.0, 0.01, write_sidarthe);

    // Days 38 - 250
    epsilon = 0.2;
    rho = 0.02;
    kappa = 0.02;
    xi = 0.02;
    sigma = 0.01;
    
    zeta = 0.025;
    eta = 0.025;
    integrate(sidarthe, x, 38.0, 350.0, 0.01, write_sidarthe); */

    return 0;
}