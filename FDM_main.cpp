/**
* Name: 		David Turner
* Description: 	Calculates the Black-Scholes Closed Form, Explicit FDM , Implicit FDM, Crank-Nicholson FDM, and Implicit SOR,
* and Crank-Nicholson SOR to price the European and American call and put options.
*
* The initial values are defined in the main function for the risk free rate, stock, sigma (volatility), dividend, Strike K, Time T,
* timestep, and step size in space. These values can be modified accordingly by the user by changing the top of the main function.
*
* ==================================================================================================================================
* To compile & run:
* 1.) Run make command in same directory as project
* 2.) Type ./FDM
*
* As follows:
* $ make
* g++ FDM_main.cpp ImplicitFDM.cpp ExplicitFDM.cpp CN_FDM.cpp \
* ImplicitSORFDM.cpp CN_SORFDM.cpp FDM_utils.cpp \
* BlackScholesFormula.cpp -o FDM
*
* $ ./FDM
*
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include "FDM_utils.h"

using namespace std;

// function prototypes defined in other files
double BlackScholesCall(double S, double K, double r, double q, double sigma, double expiry);
double BlackScholesPut(double S, double K, double r, double q, double sigma, double expiry);
double ExplicitFDM(double S, double K, double r, double q, double sigma, double expiry, double dx, double dtau, int call_or_put, int amer_or_eur);
double ImplicitFDM(double S, double K, double r, double q, double sigma, double expiry, double dx, double dtau, int call_or_put, int amer_or_eur);
double CN_FDM(double S, double K, double r, double q, double sigma, double expiry, double dx, double dtau, int call_or_put, int amer_or_eur);
double ImplicitSORFDM(double S, double K, double r, double q, double sigma, double expiry, double dx, double dtau, int call_or_put, int amer_or_eur);
double CN_SORFDM(double S, double K, double r, double q, double sigma, double expiry, double dx, double dtau, int call_or_put, int amer_or_eur);

int main() {
  double rate = 0.03;     //risk free rate
  double stock = 20.0;    //spot price of stock
  double sigma=0.8;       //volatility
  double div =0.04;       //dividend yield
  double K = 20.0;        //strike
  double T = 1.0;         //expiry
  double dtau = 0.00125;  //timestep size
  double dx = 0.05;       //step size in space

  double BSEurCall, BSEurPut;
  double ExEurCall, ExEurPut, ExAmCall, ExAmPut;
  double ImEurCall, ImEurPut, ImAmCall, ImAmPut;
  double CNEurCall, CNEurPut, CNAmCall, CNAmPut;
  double ImSOREurCall, ImSOREurPut, ImSORAmCall, ImSORAmPut;
  double CNSOREurCall, CNSOREurPut, CNSORAmCall, CNSORAmPut;

  // Values of Options are calculated from corresponding function
  BSEurCall = BlackScholesCall(stock,K,rate,div,sigma,T);
  BSEurPut = BlackScholesPut(stock,K,rate,div,sigma,T);

  ExEurCall = ExplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,1,1);
  ExEurPut = ExplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,1);
  ExAmCall = ExplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,1,0);
  ExAmPut = ExplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,0);

  ImEurCall = ImplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,1,1);
  ImEurPut = ImplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,1);
  ImAmCall = ImplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,1,0);
  ImAmPut = ImplicitFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,0);

  CNEurCall = CN_FDM(stock,K,rate,div,sigma,T,dx,dtau,1,1);
  CNEurPut = CN_FDM(stock,K,rate,div,sigma,T,dx,dtau,-1,1);
  CNAmCall = CN_FDM(stock,K,rate,div,sigma,T,dx,dtau,1,0);
  CNAmPut = CN_FDM(stock,K,rate,div,sigma,T,dx,dtau,-1,0);

  ImSOREurCall = ImplicitSORFDM(stock,K,rate,div,sigma,T,dx,dtau,1,1);
  ImSOREurPut = ImplicitSORFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,1);
  ImSORAmCall = ImplicitSORFDM(stock,K,rate,div,sigma,T,dx,dtau,1,0);
  ImSORAmPut = ImplicitSORFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,0);

  CNSOREurCall = CN_SORFDM(stock,K,rate,div,sigma,T,dx,dtau,1,1);
  CNSOREurPut = CN_SORFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,1);
  CNSORAmCall = CN_SORFDM(stock,K,rate,div,sigma,T,dx,dtau,1,0);
  CNSORAmPut = CN_SORFDM(stock,K,rate,div,sigma,T,dx,dtau,-1,0);

  // Code below is for formatting output table
  cout << endl << "Table 1: Summary of values calculated by different numeric methods" << endl << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "METHOD" << '|' << setfill(' ') << setw(14) << "OPTION TYPE" << '|' << setfill(' ') << setw(13) << "OPTION VALUE" << '|' << setfill(' ') << setw(6) << right << "ERROR" << '|' << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Close form Black Scholes" << '|' << setfill(' ') << setw(14) << "European call" << '|' << setfill(' ') << setw(13) << fixed <<  setprecision(6) << right << BSEurCall << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Close form Black Scholes" << '|' << setfill(' ') << setw(14) << "European put" << '|' << setfill(' ') << setw(13) << fixed <<  setprecision(6) << right << BSEurPut << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Explicit FDM" << '|' << setfill(' ') << setw(14) << "European call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ExEurCall << '|' << setfill(' ') << setw(5) << setprecision(2) << (ExEurCall-BSEurCall)*100/BSEurCall << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Explicit FDM" << '|' << setfill(' ') << setw(14) << "European put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ExEurPut << '|' << setfill(' ') << setw(5) << setprecision(2) << (ExEurPut-BSEurPut)*100/BSEurPut << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Explicit FDM" << '|' << setfill(' ') << setw(14) << "American call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ExAmCall << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Explicit FDM" << '|' << setfill(' ') << setw(14) << "American put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ExAmPut << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit FDM" << '|' << setfill(' ') << setw(14) << "European call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImEurCall << '|' << setfill(' ') << setw(5) << setprecision(2) << (ImEurCall-BSEurCall)*100/BSEurCall << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit FDM" << '|' << setfill(' ') << setw(14) << "European put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImEurPut << '|' << setfill(' ') << setw(5) << setprecision(2) << (ImEurPut-BSEurPut)*100/BSEurPut << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit FDM" << '|' << setfill(' ') << setw(14) << "American call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImAmCall << '|' << setfill(' ') << setw(6) << right << "N.A." << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit FDM" << '|' << setfill(' ') << setw(14) << "American put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImAmPut << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Crank-Nicholson FDM" << '|' << setfill(' ') << setw(14) << "European call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNEurCall << '|' << setfill(' ') << setw(5) << setprecision(2) << (CNEurCall-BSEurCall)*100/BSEurCall << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Crank-Nicholson FDM" << '|' << setfill(' ') << setw(14) << "European put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNEurPut << '|' << setfill(' ') << setw(5) << setprecision(2) << (CNEurPut-BSEurPut)*100/BSEurPut << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Crank-Nicholson FDM" << '|' << setfill(' ') << setw(14) << "American call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNAmCall << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Crank-Nicholson FDM" << '|' << setfill(' ') << setw(14) << "American put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNAmPut << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit (SOR) FDM" << '|' << setfill(' ') << setw(14) << "European call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImSOREurCall << '|' << setfill(' ') << setw(5) << setprecision(2) << (ImSOREurCall-BSEurCall)*100/BSEurCall << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit (SOR) FDM" << '|' << setfill(' ') << setw(14) << "European put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImSOREurPut << '|' << setfill(' ') << setw(5) << setprecision(2) << (ImSOREurPut-BSEurPut)*100/BSEurPut << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit (Projected SOR) FDM" << '|' << setfill(' ') << setw(14) << "American call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImSORAmCall << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Implicit (Projected SOR) FDM" << '|' << setfill(' ') << setw(14) << "American put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << ImSORAmPut << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Crank-Nicholson (SOR) FDM" << '|' << setfill(' ') << setw(14) << "European call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNSOREurCall << '|' << setfill(' ') << setw(5) << setprecision(2) << (CNSOREurCall-BSEurCall)*100/BSEurCall << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(35) << "Crank-Nicholson (SOR) FDM" << '|' << setfill(' ') << setw(14) << "European put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNSOREurPut << '|' << setfill(' ') << setw(5) << setprecision(2) << (CNSOREurPut-BSEurPut)*100/BSEurPut << "%" << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(24) << "Crank-Nicholson (Projected SOR) FDM" << '|' << setfill(' ') << setw(14) << "American call" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNSORAmCall << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << left << setfill(' ') << setw(24) << "Crank-Nicholson (Projected SOR) FDM" << '|' << setfill(' ') << setw(14) << "American put" << '|' << setfill(' ') << setw(13) <<  setprecision(6) << right << CNSORAmPut << '|' << setfill(' ') << setw(6) << "N.A." << '|' << endl;

  cout << '|' << right << setfill('-') << setw(36) << '|' << setfill('-') << setw(15) << '|' << setfill('-') << setw(14) << '|' << setfill('-') << setw(7) << '|' << endl;

  return 0;
};
