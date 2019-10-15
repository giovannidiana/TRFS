#include<iostream>
#include<cmath>
#include<fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

double fun(double *x){
	 return pow(x[0]-1,2)+pow(x[1]-3,2);
 }

void min_calc_rand(double (*f)(double*,double,int),double lambda,int modelID,int rand_init, double* x,double &val, const int dim,double step_init,double rate,double Tinit,string opt){
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r,rand_init); 
	
	double xnew[dim];
	double ranvec[dim];
	int i=0,j;
	int call=0;
	double inc;
	int counter=0;
	double T;
	
	ofstream outfile("par_walk.dat");
	double fMAX=40000;
	
	double f_inc,f_loc=fMAX,f_new;
	if(opt=="noinit"){
		// FIND A GOOD INITIALIZATION
		while(f_loc>=fMAX){
			for(i=0;i<dim;i++) x[i]=gsl_ran_flat(r,0,1);
			f_loc=f(x,lambda,modelID);
		}
		//for(i=0;i<dim;i++) cout<<x[i]<<' ';
		//cout<<f(x,lambda,modelID)<<endl;
	}
	
	f_loc=f(x,lambda,modelID);
    //cout<<"INIT END"<<endl;
	
	inc=10; 
	int nsteps=0;
	int time=0;
	double delta[dim];
	for(i=0;i<dim;i++) delta[i]=1;
	while(nsteps<50000){
		time++;
		T=Tinit/(time/rate);
		
		if(gsl_ran_flat(r,0,1)<0.1){
			for(i=0;i<dim;i++) xnew[i]=gsl_ran_flat(r,0,1);
		} else {
			for(i=0;i<dim;i++){
				delta[i]=min(min(step_init,x[i]),min(step_init,1-x[i]));
				//delta=step_init;
				xnew[i]=gsl_ran_flat(r,max(0.,x[i]-delta[i]),min(1.,x[i]+delta[i]));
				//if(xnew[i]<1e-100) xnew[i]=0;
			}
		}
		f_new=f(xnew,lambda,modelID);
		inc=f_new-f_loc;
		nsteps++;
		if(gsl_ran_flat(r,0,1)<exp((f_new-f_loc)/T)){
			if(inc>0) {
				//cout<<nsteps<<' '<<T<<"          \r"<<flush;
				nsteps=0;
			} 
		
		    outfile<<f_new<<' ';
			for(i=0;i<dim;i++){
				x[i]=xnew[i];
				outfile<<xnew[i]<<' ';
			}
			outfile<<endl;
			f_loc=f_new;

		} 
		
	}

	val=f_loc;
	gsl_rng_free(r);

	
}



