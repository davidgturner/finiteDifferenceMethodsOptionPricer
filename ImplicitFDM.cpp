/**
* ==================================================================================================================================
* Name: 		David Turner
* Description: Solves Black Scholes equation using Implicit finite difference method using Thomas algorithm.
*
* */

#include <cmath>
#include "FDM_utils.h"

using namespace std;

/**
 * Solves Black Scholes equation using Implicit finite difference method
 *   Matrix division is solved using tridiagonal Thomas algorithm
 * Inputs: double S (spot price)
 *         double K (strike price)
 *         double r (risk free rate)
 *         double q (dividend rate)
 *         double sigma (volatility)
 *         double expiry (time to expiry)
 *         double dx (step size in space)
 *         double dtau (step size in time)
 *         int call_or_put (+1 for call, -1 for put)
 *         int amer_or_eur (0 for European option, 1 for American)
 * Output: double value (value of option)
 */
double ImplicitFDM(double S, double K, double r, double q, double sigma, double expiry, double dx, double dtau, int call_or_put, int amer_or_eur) {
  double *a;
  double *b;
  double *fvec;

  double w;
  double *ymat;
  int i, j;
  double *t, *x;
  double t_max, t_min, x_max, x_min;
  int N, M;

  double rp = 2*r/(sigma*sigma);
  double qp = 2*(r-q)/(sigma*sigma);
  double alpha = -0.5*(qp-1);
  double beta = -0.25*(qp-1)*(qp-1) + rp;

  // x = log(S/K) ranges from -2.5 to 2.5
  // proxy for 0 to "infinity"
  x_min = -2.5;
  x_max = 2.5;
  M = 1 + ((x_max - x_min)/dx);

  x = new double[M];

  // setup x vector with M elements at steps of dx
  for(i=0; i<M; i++) {
    x[i] = x_min + i*dx;
  }

  // tau vector ranges from 0 to 0.5*sigma^2*expiry
  t_min = 0.0;
  t_max = 0.5*(sigma*sigma)*expiry;
  N = 1 + ((t_max - t_min)/dtau);

  t = new double[N];

  // setup t vector with N elements at steps of dtau
  for(j=0; j<N; j++) {
    t[j] = t_min + j*dtau;
  }

  // FDM matrix
  ymat = new double[M*N];

  for(i=0; i<M; i++) {
    // Initial condition (at tau=0)
    ymat[i] = fmax(call_or_put*(exp(0.5*x[i]*(qp+1))-exp(0.5*x[i]*(qp-1))),0.0);
  }

  w = dtau/(dx*dx);

  a = new double[3*M]; // tridiagonal matrix
                       // main diagonal stored at a[3*i+1]
                       // sub diagonal stored at a[3*i]
                       // sup diagonal stored at a[3*i+2]

  a[0+0*3] = 0.0; //unused
  a[1+0*3] = 1.0;
  a[0+1*3] = 0.0;

  for(i=1; i<M-1; i++) {
    a[2+(i-1)*3] = - w;
    a[1+ i*3] = 1.0 + 2.0 * w;
    a[0+(i+1)*3] = - w;
  }

  a[2+(M-2)*3] = 0.0;
  a[1+(M-1)*3] = 1.0;
  a[2+(M-1)*3] = 0.0; //unused

  b = new double[M];
  fvec = new double[M];

  for(j=1; j<N; j++) {
    // Boundary condition at x=-2.5
    b[0]=(call_or_put>0)?0.0:w*exp(0.5*(qp-1)*x[0]+0.25*(qp-1)*(qp-1)*t[j]);

    for(i=1; i<M-1; i++) {
      b[i] = ymat[i+(j-1)*M]; // copy current column to b
    }
    // Boundary condition at x=2.5
    b[M-1] = (call_or_put>0)?w*exp(0.5*(qp+1)*x[M-1]+0.25*(qp+1)*(qp+1)*t[j]):0.0;

    delete [] fvec; // thomas_method allocates new fvec, so we delete it first

    fvec = thomas_method(M, a, b); // solves fvec = a \ b, to get interior points

    for(i=0; i<M; i++) {
      if(amer_or_eur==1)
	ymat[i+j*M] = fmax(fvec[i],w*call_or_put*(exp(0.5*x[i]*(qp+1))-exp(0.5*x[i]*(qp-1)))); // check for early exercise
      else
	ymat[i+j*M] = fvec[i];
    }
  }

  i = M/2; // value of option at the money (S == K)
  j = N-1; // value at tau (t=0)

  double value = ymat[i+j*M]*K*exp(alpha*x[i]+beta*t[j]);

  delete [] a;
  delete [] b;
  delete [] ymat;
  delete [] t;
  delete [] x;

  return value;
}
