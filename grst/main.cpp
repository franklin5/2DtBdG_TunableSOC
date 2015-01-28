//
//  main.cpp
//  tBdG
//
//  Created by Dong Lin on 1/29/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#include "stdcpp.h"
#include "grst.h"
#include "lgwt.h"
#include "tBdG.h"
int main(){
	// *********************************** //
	// --> grst solution from gap and number equation
	// *********************************** //
		cGRST grst; // zeeman field is up to 2.416
		grst.compute();
		/*
	// TODO: initial Zeeman field as input index. Need to look up value in grst.
	int nZeeman = 2400;
	FILE *sf_input;
	sf_input = fopen ("superfluid.dat","r"); // Zeeman, Delta, Mu, Eg
	assert (sf_input != NULL);
	double delta0, mu0, Eg;
	sPara para; // E_b, Zeeman, SOC
	para.t = 0.2; para.v = 1.2; // TODO: make sure it is the same as in grst.cpp
	for (int nsf_input = 0; nsf_input <= nZeeman; ++nsf_input) {
		fscanf(sf_input, "%lf %lf %lf %lf", &para.h, &delta0, &mu0, &Eg);
	}
	fclose (sf_input);

	complex<double> delta(delta0,0);
	sPhys phys;
	phys.mu=mu0;    phys.delta=delta;
	sGauss gauss;
	//	int    NK = 1200; 	double kc_bdg = 30;
	int    NK = 4000; 	double kc_bdg = 20;
	double *gauss_k = new double [NK];double *gauss_w_k = new double [NK];
	gauss_lgwt(NK,0,kc_bdg,gauss_k,gauss_w_k);
	gauss.N=NK;    gauss.kc=kc_bdg;
	gauss.gauss_x = gauss_k;    gauss.gauss_w = gauss_w_k;

	// *********************************** //
	// --> time dependent BdG evolution
	// *********************************** //
	// --> class ctBdG for the next time step update of delta
	// !!! (note: mu is not time-dependent)

	ctBdG tBdG(phys, para, gauss);
	tBdG.Initialize_Euabv();

	para.h = 2.4;
	tBdG.quench(para.h);

	ofstream bdg_output;
      	char filename[]="test";
	//	char filename[]="tBdG_0670_dk_4000_dt_0001_long.dat";
//	char hi_buffer [100], hf_buffer[100];
//	strcat(filename,sprintf());
//	strcat(filename,".dat");
	bdg_output.open(filename); // TODO: modify file name!!!
	assert(bdg_output.is_open());
	bdg_output.precision(16);
	cout.precision(16);
	double N0, N1;
	double dt =0.001, mu;
	for (int i = 0; i < 1; ++i) {
	  // para.h = dt*i/100;    // slow change of Zeeman field.
	  // tBdG.quench(para.h);
		tBdG.update_Delta(dt, phys.delta, N0, N1);
		tBdG.get_mu(mu);
	
		  cout << "t = "<< i*dt 
		    //<< '\t'	   << "Delta = "     << abs(phys.delta) << '\t'


			//<< "mu = " << mu << endl
			//<< "real - |Delta|cos(2mu*t)= "<< phys.delta.real() - abs(phys.delta)*cos(2*mu*i*dt)<< '\t'
			//<< "imag + |Delta|sin(2mu*t) = "<< phys.delta.imag()+ abs(phys.delta)*sin(2*mu*i*dt)
//					<< "mu^inf = " << atan(phys.delta.imag() / phys.delta.real())/(2*i*dt)
				<< "n0/n = " << N0*(2*M_PI) << '\t'
					<< "n1/n = " << N1*(2*M_PI)
					<<endl;
		
			bdg_output << i*dt << '\t'
			  //					<< abs(phys.delta) << '\t'
			  //	<< abs(para.h-abs(phys.delta)) << '\t'
			  //    << para.h << '\t'
					<< phys.delta.real() << '\t' << phys.delta.imag() << '\t'
			  //	<< atan(phys.delta.imag() / phys.delta.real())/(2*i*dt) << '\t'
			  << N0*(2*M_PI) << '\t' << N1*(2*M_PI)
					<< endl;
	}
	delete []gauss_k;
	delete []gauss_w_k;
	bdg_output.close();
*/
    return 0;
}

