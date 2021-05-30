Finite Difference Methods:

https://en.wikipedia.org/wiki/Finite_difference_method

This project calculates the Black-Scholes Closed Form, Explicit FDM , Implicit FDM, Crank-Nicholson FDM, and Implicit SOR, * and Crank-Nicholson SOR to price the European and American call and put options.
*
* The initial values are defined in the main function for the risk free rate, stock, sigma (volatility), dividend, Strike K, Time T, timestep, and step size in space. These values can be modified accordingly by the user by changing the top of the main function.

=========================================================

Compile & Run Instructions / Screen Output:

To compile & run:
1.) Run make command in same directory as project
2.) Type ./FDM
As follows:

1.)
$ make
g++ FDM_main.cpp ImplicitFDM.cpp ExplicitFDM.cpp CN_FDM.cpp \
ImplicitSORFDM.cpp CN_SORFDM.cpp FDM_utils.cpp \
BlackScholesFormula.cpp -o FDM

2.) after first step it will compile to an FDM.exe file which can be executed like this
$ ./FDM

Installation / Troubleshooting Tips:
Make sure g++ and make are on the os path variable. 

For example, if you have a MingGW installation, then put this installation's bin folder on the path. 
Then make sure the make is available. 

For example:
MingGW installation: C:\DevTools\MinGW
make - mingw32-make.exe
g++ - g++

=========================================================

Author Name: David Turner

Summary and Analysis

The Black-Scholes equation is to be solved numerically and compared to the closed form equations.
For the given input parameters, i.e.
• Asset price S = 20
• Strike price K = 20
• Risk free rate r = 0.03
• Dividend yield q = 0.04
• Volatility σ = 0.8
• Maturity T = 1
The closed form equation gives values of $5.91 and $6.10 respectively for European call and put options. To
solve numerically, a grid must be defined over space and time. The space and time dimensions can be
transformed from S, t of the standard Black-Scholes to x, τ where x = log(S/K) and τ = 0.5 σ2

(T-t). The space
and time steps can then be chosen as dx = 0.05 and d τ = 0.00125. Using these values, the stability factor λ = d τ/
(dx)2
is equal to 0.5, and hence would be stable for the explicit finite difference scheme.

The equation is solved numerically using several different methods including Explicit FDM (Thomas), Implicit

FDM (Thomas), Crank-Nicholson (Thomas), Implicit FDM (SOR), Implicit FDM (Projected SOR), Crank-
Nicholson (SOR), and Crank-Nicholson (Projected SOR).

The results are summarized in the following table on the next page:

Table 1. Summary of different numerical simulation methods of solving the Black-Scholes equation

The Implicit FDM and Crank-Nicolson FDM both require the solution of a matrix equation to solve the
backward Euler step. Since the matrix is tridiagonal, this matrix equation can be solved more efficiently than for
standard matrices. For this exercise, the matrix equation is solved using Thomas algorithm and Successive Over
Relaxation (SOR), or Projected SOR for the case of American options.
The results obtained for the different option types are consistent with each other, even though there is some error
from the closed form result. The error could possibly come from the discretization step, since the step size is in
terms of x = log(S/K) rather than S. While the grid is uniform in x, the step size in S becomes larger as the ratio
S/K increases.

Sample Output (should look something like this):

PS C:\Projects\finiteDifferenceMethodsOptionPricer> mingw32-make
g++ FDM_main.cpp ImplicitFDM.cpp ExplicitFDM.cpp CN_FDM.cpp \
ImplicitSORFDM.cpp CN_SORFDM.cpp FDM_utils.cpp \
BlackScholesFormula.cpp -o FDM
PS C:\Projects\finiteDifferenceMethodsOptionPricer> .\FDM.exe   

Table 1: Summary of values calculated by different numeric methods       

|-----------------------------------|--------------|-------------|------|
|METHOD                             |OPTION TYPE   |OPTION VALUE | ERROR|
|-----------------------------------|--------------|-------------|------|
|Close form Black Scholes           |European call |     5.907000|  N.A.|
|Close form Black Scholes           |European put  |     6.100122|  N.A.|
|-----------------------------------|--------------|-------------|------|
|Explicit FDM                       |European call |     6.220736| 5.31%|
|Explicit FDM                       |European put  |     6.421977| 5.28%|
|Explicit FDM                       |American call |     6.220736|  N.A.|
|Explicit FDM                       |American put  |     6.421977|  N.A.|
|-----------------------------------|--------------|-------------|------|
|Implicit FDM                       |European call |     6.190003| 4.79%|
|Implicit FDM                       |European put  |     6.388773| 4.73%|
|Implicit FDM                       |American call |     6.156959|  N.A.|
|Implicit FDM                       |American put  |     6.352955|  N.A.|
|-----------------------------------|--------------|-------------|------|
|Crank-Nicholson FDM                |European call |     6.182996| 4.67%|
|Crank-Nicholson FDM                |European put  |     6.380880| 4.60%|
|Crank-Nicholson FDM                |American call |     6.182325|  N.A.|
|Crank-Nicholson FDM                |American put  |     6.380138|  N.A.|
|-----------------------------------|--------------|-------------|------|
|Implicit (SOR) FDM                 |European call |     6.190003| 4.79%|
|Implicit (SOR) FDM                 |European put  |     6.388753| 4.73%|
|Implicit (Projected SOR) FDM       |American call |     6.156958|  N.A.|
|Implicit (Projected SOR) FDM       |American put  |     6.352934|  N.A.|
|-----------------------------------|--------------|-------------|------|
|Crank-Nicholson (SOR) FDM          |European call |     6.179574| 4.61%|
|Crank-Nicholson (SOR) FDM          |European put  |     6.377203| 4.54%|
|Crank-Nicholson (Projected SOR) FDM|American call |     6.164172|  N.A.|
|Crank-Nicholson (Projected SOR) FDM|American put  |     6.360504|  N.A.|
|-----------------------------------|--------------|-------------|------|