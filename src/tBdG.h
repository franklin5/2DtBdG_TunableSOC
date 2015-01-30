//
//  tBdG.h
//  tBdG
//
//  Created by Dong Lin on 1/30/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#ifndef tBdG_tBdG_h
#define tBdG_tBdG_h

#include "stdcpp.h"
#include "dist.h"
class ctBdG: public cDistribute{
private:
  double _Eb, _hi, _omega, _v, _mu, _dt, _total_t, _Omega1;
    complex<double> _delta;
    double _N0, _N1;
    int _NK, _NK2;
    double _kc, _Ueff;
    double *_gauss_k, *_gauss_w_k;
    complex<double> myI;
    MatrixXd _bdg_E;
    MatrixXcd _bdg, _bdg_u, _bdg_a, _bdg_b, _bdg_v;
    double* local_Delta_k_r, *local_Delta_k_i, *total_Delta_k_r, *total_Delta_k_i;
public:
    ctBdG (const int rank, const int size, const int root) : cDistribute(rank,size,root){}
    ~ctBdG(){
    	delete []_gauss_k;
    	delete []_gauss_w_k;
	/*    	delete []local_Delta_k_r;
    	delete []local_Delta_k_i;
    	if (_rank==_root) {
    		delete []total_Delta_k_r;
    		delete []total_Delta_k_i;
		}*/
    }
    void input();
    complex<double> DELTA_K(int,int);
    void RK_Propagator(int, int);
    void compute_DeltaK(complex<double>&, double&, double&, double);
    double N0_CF(int, int);
    double N1_CF(int, int);
    void tuning();
    void Initialize_Euabv();
    void construct_BdG(double, double, double);
    double spintexture(double);
};





#endif

