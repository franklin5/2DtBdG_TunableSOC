/*
 * grst.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ld7
 */
#include "grst.h"
void cGRST :: compute(){
  double omega = 0.1, hi = 1.0, Omega1 = 0.6;
  char filename[150];
  sprintf(filename,"superfluid_instant_hi_%gOmega1_%g_omega_%gomega.dat",hi,Omega1,omega);
	ofstream superfluid_output;
	superfluid_output.open(filename);
	superfluid_output.is_open();

	int    NK = 500;    double kc = 10;
	double *gauss_k = new double [NK];double *gauss_w_k = new double [NK];
	gauss_lgwt(NK,0,kc,gauss_k,gauss_w_k);
	sGauss gauss;
	gauss.N=NK; gauss.kc=kc;
	gauss.gauss_x = gauss_k;gauss.gauss_w = gauss_w_k;
	double mu0= 0.05;    double delta0 = 0.05; double Eg;
	double Zeeman = 1.3, soc = 1.2;
	int nt = 0;double dt = 0.1;
	do {
	
	  // *********************************** //
	  // --> input parameters of the system.
	  sPara para = {0.2, Zeeman, soc};
	  // *********************************** //
	  cSeek_Gap_Number GapNumber(delta0, mu0 ,para,gauss);
	  GapNumber.gapnumber();
	  GapNumber.getresult(delta0, mu0, Eg);
	  cout << nt*dt << endl;
	  /*	  cout.precision(16);
	  cout << "Zeeman = " << Zeeman << endl;
	  cout << "Delta = " << delta0 << endl;
	  cout << "Mu = "    << mu0    << endl;*/
	  GapNumber.printEg();
	  superfluid_output.precision(16);
	  superfluid_output << Zeeman << '\t' << delta0 << '\t'
			    << mu0 << '\t' << Eg << '\t' << nt*dt << endl;
	  nt += 1;
	  Zeeman = hi+Omega1/2*cos(omega*nt*dt);
	} while(nt*dt < 100);
	delete []gauss_k;
	delete []gauss_w_k;
	superfluid_output.close();
}


