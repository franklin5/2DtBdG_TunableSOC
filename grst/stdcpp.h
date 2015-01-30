//
//  stdcpp.h
//  tBdG
//
//  Created by Dong Lin on 1/29/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#ifndef tBdG_stdcpp_h
#define tBdG_stdcpp_h

#include <iostream>
#include <new>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <Eigen/Eigenvalues>
using namespace std;
using namespace Eigen;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct sPhys {
    double mu;
    complex<double> delta;
};

struct sPara {
    double t; // inverse scattering length, or bound state energy Eb in 2D
    double h; // zeeman field
    double v; // rashba SOC strength
};

struct sGauss {
    int N;
    double kc;
    double *gauss_x;
    double *gauss_w;
};


#endif
