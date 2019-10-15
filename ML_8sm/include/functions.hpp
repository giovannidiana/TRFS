#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<functional>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "constants.hpp"
#include<boost/tuple/tuple.hpp>

//#include "model_XTF.hpp"

using namespace std;
using namespace constants;

typedef boost::tuple<int,double> idouble;

class record {
	public:
	int parent;
	double t[NPAR];
	int model;
	double chi2;
	double w;
	int alpha;
};


struct SortBy{
        inline bool operator() (const idouble& t1, const idouble& t2){
	            return (t1.get<1>() > t2.get<1>());
        }
};


void normalize(int size, double* P){
	double norm=0;
	for(int i=0; i<size;i++) norm+=P[i];
	for(int i=0; i<size;i++) P[i]/=norm;
}

template <typename MODEL, typename STATE>
void single_step(double &t, STATE &s, MODEL &model, int &next_proc, gsl_rng* r){
	const int NPROC=model.NPROC;
	double P[NPROC];
	model.SetProb(P,s,r,0);
	double restime=0;
	for(int i=0;i<NPROC;i++) restime+=P[i];
	normalize(NPROC,P);
	
	t+=gsl_ran_exponential(r,1./restime);
	gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(NPROC, P);
	next_proc=gsl_ran_discrete(r,g);
	model.move(next_proc,s);
	gsl_ran_discrete_free(g);
}

int get_ID(double* params){
	int bin[9];
	int id=0;
	int k;
	bin[0]=params[0] > 1 ? 0 : 1;
	bin[1]=params[1] > 1 ? 0 : 1;
	bin[2]=params[2] > 1 ? 0 : 1;
	bin[3]=params[3] > 1 ? 0 : 1;
	bin[4]=params[4] > 1 ? 0 : 1;
	bin[5]=params[5] > 1 ? 0 : 1;
	bin[6]=params[6] > 1 ? 0 : 1;
	bin[7]=params[7] > 1 ? 0 : 1;
	bin[8]=params[8] > 1 ? 0 : 1;

    for(k=0;k<9;k++) id+=pow(2,k)*bin[k];

	return(id);
}

void get_arch(int ID, int* arch){
	int num=ID;
	arch[0]=(num%2==0) ? 0 : 1;
	for(int i=1;i<4;i++){
		num=(num-arch[i-1])/2;
		arch[i]=(num%2==0) ? 0 : 1;
	}
}

void get_top(double* Lvec, int* particle_model, int (&top)[NMOD][NTOP], int* load, double (&weights_top)[NMOD][NTOP], int w_size){
	int i,k;
	idouble * weights_sort = new idouble[w_size];
	int model_ID, particle_ID;

	for(i=0;i<w_size;i++){
		weights_sort[i].get<0>()=i;
		weights_sort[i].get<1>()=Lvec[i];
	}

	for(i=0;i<NMOD;i++) load[i]=0;

	sort(weights_sort,weights_sort+w_size,SortBy());

	for(k=0;k<w_size;k++){
		particle_ID=weights_sort[k].get<0>();
		model_ID=particle_model[particle_ID];
		if(load[model_ID]<NTOP){
			top[model_ID][load[model_ID]]=particle_ID;
			weights_top[model_ID][load[model_ID]]=weights_sort[k].get<1>();
			load[model_ID]++;
		}
	}

	for(k=0;k<NMOD;k++){
		if(load[k]>0) normalize(load[k],weights_top[k]);
	}
		
}
				

		
#endif
