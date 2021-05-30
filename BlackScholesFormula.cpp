/**
* ==================================================================================================================================
* Name: 		David Turner
* Description: calculates the closed form solution for black scholes call/put option price.
* */

#include <cmath>

using namespace std;

/*----------------------------------------------------------------------------------------
    Normal Density Function
 -----------------------------------------------------------------------------------------*/
double NDensity(double x)
{
    double f;
    double Pi = 4.0*atan(1.0);
    f = (1.0/sqrt(2.0*Pi)*exp(-x*x/2.0));
    return f;
};

// standard normal cumulative distribution function
double CNDist(double x) {
    double a1 = 0.319381530;
    double a2 = -0.356563782;
    double a3 = 1.781477937;
    double a4 = -1.821255978;
    double a5 = 1.330274429;

    double z = 1.0/(1.0 + 0.2316419*abs(x));
    double right_area = NDensity(x)*z*((((a5*z + a4)*z + a3)*z + a2)*z + a1);

    double F;

    if(x>=0)
    {
        F = 1.0 - right_area;
    }
    else
    {
        F = right_area;
    }
    return F;
};

/**
*   Analytical solution to European Call
*/
double BlackScholesCall(double S, double K, double r, double q, double sigma, double expiry) {
	double d1 = (log(S/K)+(r-q+sigma*sigma/2.0)*(expiry))/(sigma*sqrt(expiry));
	double d2 = d1 - sigma*sqrt(expiry);
    return S*exp(-q*(expiry))*CNDist(d1) - K*exp(-r*(expiry))*CNDist(d2);
};

/**
*    Analytical solutiion to European Put
*/
double BlackScholesPut(double S, double K, double r, double q, double sigma, double expiry) {
  	double d1 = (log(S/K)+(r-q+sigma*sigma/2.0)*(expiry))/(sigma*sqrt(expiry));
  	double d2 = d1 - sigma*sqrt(expiry);
  	return -S*exp(-q*(expiry))*CNDist(-d1) + K*exp(-r*(expiry))*CNDist(-d2);
};
