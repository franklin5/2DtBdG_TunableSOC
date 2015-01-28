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

class ctBdG {
private:
    double _Eb, _h, _v;
    complex<double> _delta;
    double _mu;
    int _NK;
    double _kc, _Ueff;
    double *_gauss_k, *_gauss_w_k;
    MatrixXcd _bdg_E, _bdg_u, _bdg_a, _bdg_b, _bdg_v;
public:
    ctBdG(const sPhys& phys, const sPara& para, const sGauss& gauss)
    :_Eb(para.t), _h(para.h), _v(para.v),
    _mu(phys.mu),
    _NK(gauss.N), _kc(gauss.kc),
    _Ueff ((-8.0*M_PI)/log(1.0+2.0*_kc*_kc/_Eb)) ,
    _gauss_k(gauss.gauss_x), _gauss_w_k(gauss.gauss_w),
    _bdg_E(_NK, 4), _bdg_u(_NK,4), _bdg_a(_NK,4),
    _bdg_b(_NK,4),_bdg_v(_NK,4){
  	_delta = phys.delta;
    } // this type of initialization requires ordering, same as private member orders.
    void quench(double h) {_h = h;}
    void Initialize_Euabv();
    ~ctBdG(){}
    void construct_BdG(MatrixXcd& , double, double);
    void update_Delta(double, complex<double>&, double&, double&);
    complex<double> integrand(int);
    double N0_CF(int);
    double N1_CF(int);
    void RK_Propagator(double);
    void get_mu(double& mu){mu = _mu;}
};





#endif

