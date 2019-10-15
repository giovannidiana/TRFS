#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "functions.hpp"
#include "constants.hpp"
#include<armadillo>

using namespace std;
//using namespace constants;

// Definition of the model used.
//typedef model_XTF MODEL;

// ---------------------------------------------------------------------
// As long as the model class is well defined this is model independent.
void solveCTM(double* params,int KOA,int KOB,double* out_array){

    arma::vec prob(8);
	arma::vec b(8);
	arma::mat W(8,8);

	double wxx=params[0];
	double wxy=params[1];
	double wxz=params[2];
	double wyx=params[3];
	double wyy=params[4];
	double wyz=params[5];
	double wzx=params[6];
	double wzy=params[7];
	double wzz=params[8];

	double kon=1;

	if(KOA==1){ 
		wxx=1; wxy=1; wxz=1;
		wzx=1; wzy=1; wzz=1;
	}
	if(KOB==1){ wyx=1;wyy=1;wyz=1;}
	
	W<< -3*kon << wxx  << wyy  << 0    << wzz  << 0     << 0     << 0<<arma::endr
	 <<    kon << -2*kon-wxx << 0  << wxy*wyy   << 0  << wxz*wzz   << 0     << 0<<arma::endr
	 <<    kon << 0  << -wyy-2*kon << wyx*wxx   << 0  << 0     << wyz*wzz   << 0<<arma::endr
	 <<   0  << kon  << kon  << -kon-wxx*wyx-wyy*wxy << 0  << 0     << 0     << wxz*wyz*wzz<<arma::endr
	 <<    kon << 0  << 0  << 0    << -wzz-2*kon << wzx*wxx   << wzy*wyy   << 0<<arma::endr
	 <<   0  << kon  << 0  << 0    << kon  << -kon-wzx*wxx-wxz*wzz << 0     << wyy*wxy*wzy<<arma::endr
	 <<   0  << 0  << kon  << 0    << kon  << 0     << -wzz*wyz-wyy*wzy-kon << wxx*wyx*wzx<<arma::endr
	 <<   1  << 1  << 1  << 1    << 1  << 1    << 1    << 1 <<arma::endr;

	b<<0<<0<<0<<0<<0<<0<<0<<1;

	prob=arma::solve(W,b);
	
	out_array[0]=prob(1)+prob(3)+prob(5)+prob(7);
	out_array[1]=prob(3)+prob(2)+prob(7)+prob(6);
	out_array[2]=prob(4)+prob(5)+prob(6)+prob(7);
}
// ---------------------------------------------------------------------

