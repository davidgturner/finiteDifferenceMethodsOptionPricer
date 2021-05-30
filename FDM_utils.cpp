/**
* ==================================================================================================================================
* Name: 		David Turner
* Description: a utility class encapsulating thomas method & sor method.
*
* */

#include "FDM_utils.h"
#include <cmath>

using namespace std;

/**
 * Function to implement Thomas algorithm for tridiagonal matrix
 * Inputs : int n (number of steps)
 *          double* a (length=3*n), stores tridiagonal matrix
 *          double* b (length=n), stores right hand side
 * Output : double* x (length=n), where a*x = b
 *          Thus, x = a \ b
 */
double* thomas_method(int n, double *a, double *b) {
  int i;
  double *x, *diag_new, *b_new;

  x = new double[n];
  diag_new = new double[n];
  b_new = new double[n];

  // main diagonal is stored in a[3*i+1]
  // sub diagonal is stored in a[3*i]
  // sup diagonal is stored in a[3*i+2]
  diag_new[0] = a[3*0+1];
  b_new[0] = b[0];
  // Thomas algorithm forward iteration
  for(i=1; i<n; i++) {
    diag_new[i] = a[3*i+1] - a[3*(i-1)+2]*a[3*i]/diag_new[i-1];
    b_new[i] = b[i] - b_new[i-1]*a[3*i]/diag_new[i-1];
  }

  x[n-1] = b_new[n-1]/a[3*(n-1)+1];
  // Thomas algorithm backward iteration
  for(i=n-2; i>=0; i--) {
    x[i] = (b_new[i] - a[3*i+2]*x[i+1])/diag_new[i];
  }
  delete [] diag_new;
  delete [] b_new;
  return x;
}

/**
 * Function to implement Successive OverRelaxation for tridiagonal matrix
 * Inputs : int n (number of steps)
 *          double* a (length=3*n), stores tridiagonal matrix
 *          double* b (length=n), stores right hand side
 *          double relax, the relaxation parameter
 *          int max_iter, maximum iterations before SOR terminates
 * Output : double* x (length=n), where a*x = b
 *          Thus, x = a \ b
 */
double* sor_method(int n, double *a, double *b, double relax, int max_iter) {
  int iter, i;
  double *x, *x_new;
  double square_sum = 0.0, tol = 1e-6;

  x = new double[n];
  x_new = new double[n];

  for(i=0; i<n; i++) {
    x[i] = 0.0; // initialize x vector to 0
  }

  // main diagonal is stored in a[3*i+1]
  // sub diagonal is stored in a[3*i]
  // sup diagonal is stored in a[3*i+2]
  for(iter=1;iter<=max_iter;iter++) {
    for(i=1; i<n-1; i++) {
      x_new[i] = b[i];
      // since matrix is tridiagonal, only need to update 3 elements
      x_new[i] = x_new[i] - a[3*i] * x_new[i-1];
      x_new[i] = x_new[i] - a[3*i+2] * x[i+1];
      x_new[i] = x_new[i] / a[3*i+1];
      x_new[i] = ( 1.0 - relax ) * x[i] + relax * x_new[i];
      // calculate vector norm of error
      square_sum += (x_new[i]-x[i])*(x_new[i]-x[i]);
    }
    // if change from x to x_new is below tolerance, or if
    // max iterations is reached, terminate and return best guess
    if (sqrt(square_sum) < tol || iter==max_iter) {
      delete [] x_new;
      return x;
    }
    for(i=1; i<n-1; i++) {
      x[i] = x_new[i]; // set x to x_new and continue iteration
    }
    square_sum = 0.0; // reset error back to 0 before next loop
  }

}
