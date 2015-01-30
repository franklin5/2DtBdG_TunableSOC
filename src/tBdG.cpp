/*
 * tBdG.cpp
 *
 *  Created on: Nov 11, 2014
 *      Author: ld7
 */

#include "tBdG.h"
#include "lgwt.h"

void ctBdG :: input(){
  for (int ig = 0; ig < _size; ++ig) {
    if (ig ==_rank){
      char dummyname[100];
      double dummyvalue;
      int intdummyvalue;
      FILE *input;
      input = fopen("input.txt","r");
      assert(input != NULL);
      if (ig == _root)  cout << "Starting to read in parameters from file input.txt" << endl;
      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _hi = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _hi << endl;
      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _Omega1 = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _Omega1 << endl;

      FILE *sf_input;
      sf_input = fopen ("superfluid.dat","r"); // Zeeman, Delta, Mu, Eg
      double dummyH, dummyD, dummyM, dummyE;
      assert (sf_input != NULL);
      for (int nsf_input = 0; nsf_input <= int((_hi+_Omega1/2)/0.001); ++nsf_input) {
	fscanf(sf_input, "%lf %lf %lf %lf", &dummyH, &dummyD, &dummyM, &dummyE);
      }
      fclose (sf_input);
      _delta.real() = dummyD; _delta.imag() = 0.0;
      _N0 = 0.0; _N1 = 0.0;
      _mu = dummyM;

      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _omega = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _omega << endl;
      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _Eb = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _Eb << endl;
      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _v = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _v << endl;
      fscanf(input,"%s %d", dummyname, &intdummyvalue);
      _NK = intdummyvalue;    if (ig == _root) cout << dummyname << "=" << _NK << endl;
      _NK2 = _NK;
      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _kc = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _kc << endl;
      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _dt = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _dt << endl;
      fscanf(input,"%s %lf", dummyname, &dummyvalue);
      _total_t = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _total_t << endl;


      fclose(input);

      _Ueff = (-8.0*M_PI)/log(1.0+2.0*_kc*_kc/_Eb);

      _bdg.resize(4,4);
      _gauss_k = new double [_NK];
      _gauss_w_k = new double [_NK];
      gauss_lgwt(_NK,0,_kc,_gauss_k,_gauss_w_k);

      myI = complex<double> (0.0,1.0);
    }
  }
}

void ctBdG:: Initialize_Euabv(){
	SelfAdjointEigenSolver<MatrixXcd> ces;
	distribution(_NK2);
    int nk;
    for (int ig = 0; ig < _size; ++ig) {
		if (ig==_rank) {
			_bdg_E.resize(recvcount,4);
			_bdg_u.resize(recvcount,4);
			_bdg_a.resize(recvcount,4);
			_bdg_b.resize(recvcount,4);
			_bdg_v.resize(recvcount,4);
//			local_Delta_k_r = new double[recvcount];
//			local_Delta_k_i = new double[recvcount];
			for (int i = 0; i < recvcount; ++i) {

				nk = recvbuf[i];

				construct_BdG(_gauss_k[nk], _mu, _hi);
				ces.compute(_bdg);

				_bdg_E.row(i) = ces.eigenvalues();
				_bdg_u.row(i) = ces.eigenvectors().row(0);
				_bdg_a.row(i) = ces.eigenvectors().row(1);
				_bdg_b.row(i) = ces.eigenvectors().row(2);
				_bdg_v.row(i) = ces.eigenvectors().row(3);
			}
		}
	}
}

void ctBdG:: construct_BdG(double kx, double mu, double h){
	double xi = kx*kx-mu;
	_bdg(0,0) = complex<double> (xi+h,0.0);
	_bdg(0,1) = complex<double> (_v*kx,0.0);
	_bdg(0,2) = complex<double> (0.0,0.0);
	_bdg(0,3) = -_delta;
	_bdg(1,0) = complex<double> (_v*kx,0.0);
	_bdg(1,1) = complex<double> (xi-h,0.0);
	_bdg(1,2) = _delta;
	_bdg(1,3) = complex<double> (0.0,0.0);
	_bdg(2,0) = complex<double> (0.0,0.0);
	_bdg(2,1) = conj(_delta);
	_bdg(2,2) = complex<double> (-(xi+h),0.0);
	_bdg(2,3) = complex<double> (_v*kx,0.0);
	_bdg(3,0) = -conj(_delta);
	_bdg(3,1) = complex<double> (0.0,0.0);
	_bdg(3,2) = complex<double> (_v*kx,0.0);
	_bdg(3,3) = complex<double> (-xi+h,0.0);

}

void ctBdG::tuning(){
// tuning takes place by changing h from h_i to h_i\times\cos(omega \times t) value at all times t = 0.
	complex<double> localDelta;
	double localN0, localN1;
	ofstream delta_output;
//	ofstream akx;
//	ofstream aky, Delta_K_r, Delta_K_i;
	if (_rank == _root) {
//		total_Delta_k_r = new double[_NK2];
//		total_Delta_k_i = new double[_NK2];
/*		string filename = "hi_" + to_string(_hi) + "hf_" + to_string(_hf) + ".dat"; // --> This is only available for c++11 feature...
		delta_output.open(filename.c_str());*/
		char filename[50];
		sprintf(filename,"hi_%gOmega1_%gomega_%g.dat",_hi,_Omega1,_omega); // Use the shortest representation %g
		delta_output.open(filename);delta_output.is_open();
//		akx.open("akx.OUT");akx.is_open();
//		aky.open("aky.OUT");aky.is_open();
		delta_output.precision(16);
//		akx.precision(16);
//		aky.precision(16);
//		for (int nk = 0; nk < _NK; ++nk) {
//			akx << _gauss_k[nk] << endl;
//			aky << _gauss_k[nk] << endl;
//		}
//		akx.close();
//		aky.close();
//		sprintf(filename,"hi_%gomega_%g_Delta_K_r.dat",_hi,_omega);
//		Delta_K_r.open(filename);Delta_K_r.is_open();
//		sprintf(filename,"hi_%gomega_%g_Delta_K_i.dat",_hi,_omega);
//		Delta_K_i.open(filename);Delta_K_i.is_open();
//		cout.precision(16);
//		cout << _delta << endl;
		//		delta_output << 0.0 << '\t' << _delta.real() << '\t' << _delta.imag() << '\t' << _hi << endl;
//		Delta_K_r.precision(16);Delta_K_i.precision(16);
	}
	double ht, sz0;
	for (int nt = 1; nt < int(_total_t/_dt); ++nt) {
		ht = _hi+_Omega1/2*cos(_omega*nt*_dt);
		compute_DeltaK(localDelta, localN0, localN1, ht);
		// gather or reduce to rank 0 to update \Delta(t) value
//		MPI_Reduce(&localDelta, &_delta, 1, MPI_C_DOUBLE_COMPLEX, MPI_SUM, _root, MPI_COMM_WORLD);
		MPI_Reduce(&localDelta.real(), &_delta.real(), 1, MPI_DOUBLE, MPI_SUM, _root, MPI_COMM_WORLD);
		MPI_Reduce(&localDelta.imag(), &_delta.imag(), 1, MPI_DOUBLE, MPI_SUM, _root, MPI_COMM_WORLD);
		// scatter or broadcast \Delta(t) to all rank
//		MPI_Bcast(&_delta, 1, MPI_C_DOUBLE_COMPLEX, _root, MPI_COMM_WORLD);
		MPI_Bcast(&_delta.real(), 1, MPI_DOUBLE, _root, MPI_COMM_WORLD);
		MPI_Bcast(&_delta.imag(), 1, MPI_DOUBLE, _root, MPI_COMM_WORLD);
		MPI_Reduce(&localN0, &_N0, 1, MPI_DOUBLE, MPI_SUM, _root, MPI_COMM_WORLD);
		MPI_Reduce(&localN1, &_N1, 1, MPI_DOUBLE, MPI_SUM, _root, MPI_COMM_WORLD);
		// Gather \Delta(k,t) complex function
//		MPI_Gatherv(local_Delta_k_r,recvcount, MPI_DOUBLE, total_Delta_k_r,recvcounts,displs_r,MPI_DOUBLE, _root, COMM_WORLD);
//		MPI_Gatherv(local_Delta_k_i,recvcount, MPI_DOUBLE, total_Delta_k_i,recvcounts,displs_r,MPI_DOUBLE, _root, COMM_WORLD);
		// rinse and repeat
		if (_rank == _root) {
		  sz0 = spintexture(ht);
//			if (nt%int(0.1/_dt)==0) {
		  if (nt%1000==0) {
		    cout << nt*_dt << '\t' << _delta << '\t' << sz0 << endl;
		  }
		  delta_output << nt*_dt << '\t'
			       << _delta.real() << '\t' << _delta.imag() << '\t'
			       << ht << '\t' << 2*_N0*(2*M_PI) << '\t' << 2*_N1*(2*M_PI) << '\t' << sz0 << endl;
		  //				 for (int nk = 0; nk < _NK2; ++nk) {
		  //					Delta_K_r << total_Delta_k_r[nk] << '\t';
		  //					Delta_K_i << total_Delta_k_i[nk] << '\t';
		  //				}
		  //				Delta_K_r << endl;Delta_K_i << endl;
//			}
		}
	}
	//	Delta_K_r.close();Delta_K_i.close();
	delta_output.close();
}

double  ctBdG::spintexture(double ht){
  double kx = 0.0, mu = 0.0, xi = kx*kx-mu;
  Matrix4cd _bdg_H;
  _bdg_H(0,0) = complex<double> (xi+ht,0.0);
  _bdg_H(0,1) = complex<double> (_v*kx,0.0);
  _bdg_H(0,2) = complex<double> (0.0,0.0);
  _bdg_H(0,3) = -_delta;
  _bdg_H(1,0) = complex<double> (_v*kx,0.0);
  _bdg_H(1,1) = complex<double> (xi-ht,0.0);
  _bdg_H(1,2) = _delta;
  _bdg_H(1,3) = complex<double> (0.0,0.0);
  _bdg_H(2,0) = complex<double> (0.0,0.0);
  _bdg_H(2,1) = conj(_delta); // should have been conjugate of Delta here. Note for the tBdG update!                            
  _bdg_H(2,2) = complex<double> (-(xi+ht),0.0);
  _bdg_H(2,3) = complex<double> (_v*kx,0.0);
  _bdg_H(3,0) = -conj(_delta); // should have been conjugate of Delta here. Note for the tBdG update!                                                           
  _bdg_H(3,1) = complex<double> (0.0,0.0);
  _bdg_H(3,2) = complex<double> (_v*kx,0.0);
  _bdg_H(3,3) = complex<double> (-xi+ht,0.0);
  
  SelfAdjointEigenSolver<MatrixXcd> ces;
  ces.compute(_bdg_H);
  double tmpSz0 = pow(abs(ces.eigenvectors().row(2)[2]),2.0) - pow(abs(ces.eigenvectors().row(3)[2]),2.0);
  tmpSz0 += pow(abs(ces.eigenvectors().row(2)[3]),2.0) - pow(abs(ces.eigenvectors().row(3)[3]),2.0);
  //  cout << tmpSz0 << endl;
  return tmpSz0; 
}


void ctBdG:: compute_DeltaK(complex<double>& localDelta, double& localN0, double& localN1, double ht){
    int nk;
    complex<double> result;
    double N0_tmp, N1_tmp;
	for (int ig = 0; ig < _size; ++ig) {
		if (ig==_rank) {
			localDelta = complex<double> (0.0,0.0);
			localN0 = 0; localN1 = 0;
			for (int i = 0; i < recvcount; ++i) {
				nk = recvbuf[i];
				construct_BdG(_gauss_k[nk], _mu, ht); // set mu = 0
				result = complex<double> (0.0,0.0);
				N0_tmp = 0.0;N1_tmp = 0.0;
				for (int eta = 0; eta < 4; ++eta) {
					RK_Propagator(i,eta);
					result += DELTA_K(i,eta);
					N0_tmp += N0_CF(i,eta);
					N1_tmp += N1_CF(i,eta);
				}
//				local_Delta_k_r[i] = result.real();
//				local_Delta_k_i[i] = result.imag();
				localDelta += -_Ueff/(4.0*M_PI)*_gauss_k[nk]*_gauss_w_k[nk]*result;
				localN0 += 1.0/(8.0*M_PI)*_gauss_w_k[nk] *_gauss_k[nk]*N0_tmp;
				localN1 += 1.0/(8.0*M_PI)*_gauss_w_k[nk] *_gauss_k[nk]*N1_tmp;
			}
			//cout << "rank " << ig << " has localDelta = " << localDelta << endl;
		}
	}
}

complex<double> ctBdG::DELTA_K(int i,int eta){
	complex<double> result;
	switch (sgn(_bdg_E(i,eta))) {
		case -1:
			result =  _bdg_u(i,eta) * conj(_bdg_v(i,eta));
			break;
		case 1:
			result =  _bdg_a(i,eta) * conj(_bdg_b(i,eta));
			break;
		case 0:
			result =  (_bdg_u(i,eta) * conj(_bdg_v(i,eta)) +
					_bdg_a(i,eta) * conj(_bdg_b(i,eta)))/2.0;
			break;
		default:
			break;
	}
	return result;
}

double ctBdG::N0_CF(int i,int eta){
  complex<double> result (0.0,0.0);
  complex<double> u, a, b, v, E;

    E = _bdg_E(i,eta);
    u = _bdg_u(i,eta);
    a = _bdg_a(i,eta);
    b = _bdg_b(i,eta);
    v = _bdg_v(i,eta);

    switch (sgn(E.real())) {
    case -1:
      result += u * conj(v);
      break;
    case 1:
      result += a * conj(b);
      break;
    case 0:
      result += (u * conj(v) + a * conj(b))/2.0;
      break;
    default:
      break;
    }
  return  pow(abs(result),2.0);
}

double ctBdG::N1_CF(int i,int eta){
  complex<double> result (0.0,0.0);
  complex<double> u, a, b, v, E;
  E = _bdg_E(i,eta);
  u = _bdg_u(i,eta);
  a = _bdg_a(i,eta);
  b = _bdg_b(i,eta);
  v = _bdg_v(i,eta);
  switch (sgn(E.real())) {
  case -1:
    result += u * conj(b);
    break;
  case 1:
    result += a * conj(v);
    break;
  case 0:
    result += (u * conj(b) + a * conj(v))/2.0;
    break;
  default:
    break;
  }
  return  pow(abs(result),2.0);
}


void ctBdG::RK_Propagator(int i, int eta){
	VectorXcd wvVec(4), k1(4), k2(4), k3(4), k4(4);
	double normwv;
	wvVec(0) = _bdg_u(i,eta);
	wvVec(1) = _bdg_a(i,eta);
	wvVec(2) = _bdg_b(i,eta);
	wvVec(3) = _bdg_v(i,eta);
	k1 = -myI * (_bdg *  wvVec);
	k2 = -myI * (_bdg * (wvVec+_dt*k1*0.5));
	k3 = -myI * (_bdg * (wvVec+_dt*k2*0.5));
	k4 = -myI * (_bdg * (wvVec+_dt*k3));
	wvVec += _dt*(k1+2.0*k2+2.0*k3+k4)/6.0;
	normwv = wvVec.norm();
	wvVec /= normwv;
	//		cout << normwv << endl;
	_bdg_u(i,eta) = wvVec(0);
	_bdg_a(i,eta) = wvVec(1);
	_bdg_b(i,eta) = wvVec(2);
	_bdg_v(i,eta) = wvVec(3);
}
