#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstring>
#include<vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "include/functions.hpp"
#include "include/solveCTM.hpp"
#include "include/constants.hpp"

using namespace std;
using namespace constants;

// ---------------------------------------------------------------------
// As long as the model class is well defined this is model independent.
int main(int argc, char** argv){
   
	int i,j,k,l;

	// Assign sample size
	int particle_ID;
    int counter;

	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r,atoi(argv[3])); 
 
    // set parameter scales
	double params[NPAR];
	double theta[NPAR];
	double likelihood,likelihood_max=1e10;
    int model;

	double out_array_mc[3];
	double out_array_mc_KOA[3];
	double out_array_mc_KOB[3];
	double out_array_mc_KOAB[3];

	vector<double>  feat_vec(9);

	double I[3][3], I_sd[3][3];
	ostringstream fname, fname_sd;
	ostringstream out_filename;
	fname<<"intervals/interval_T"<<argv[1]<<"_"<<argv[2]<<".dat";
	fname_sd<<"intervals/sd_T"<<argv[1]<<"_"<<argv[2]<<".dat";
    out_filename<<"results/T"<<argv[1]<<"_"<<argv[2]<<".dat";

	ifstream interval_file(fname.str().c_str());
	ofstream out_file(out_filename.str().c_str());

	interval_file>>I[0][0]>>I[0][1]>>I[0][2];
	interval_file>>I[1][0]>>I[1][1]>>I[1][2];
	interval_file>>I[2][0]>>I[2][1]>>I[2][2];
	
	ifstream sd_file(fname_sd.str().c_str());
	sd_file>>I_sd[0][0]>>I_sd[0][1]>>I_sd[0][2];
	sd_file>>I_sd[1][0]>>I_sd[1][1]>>I_sd[1][2];
	sd_file>>I_sd[2][0]>>I_sd[2][1]>>I_sd[2][2];

	cout<<"Ratio data:"<<endl;
	cout<<setw(19)<<"ADF"<<setw(21)<<"ASI"<<"NSM"<<endl;
	cout<<"QL404 "<<setw(10)<<I[0][0]<<' '<<setw(10)<<I_sd[0][0]<<" "<<setw(10)<<I[1][0]<<' '<<setw(10)<<I_sd[1][0]<<" "<<setw(10)<<I[2][0]<<' '<<setw(10)<<I_sd[2][0]<<endl;
	cout<<"QL402 "<<setw(10)<<I[0][1]<<' '<<setw(10)<<I_sd[0][1]<<" "<<setw(10)<<I[1][1]<<' '<<setw(10)<<I_sd[1][1]<<" "<<setw(10)<<I[2][1]<<' '<<setw(10)<<I_sd[2][1]<<endl;
	cout<<"QL435 "<<setw(10)<<I[0][2]<<' '<<setw(10)<<I_sd[0][2]<<" "<<setw(10)<<I[1][2]<<' '<<setw(10)<<I_sd[1][2]<<" "<<setw(10)<<I[2][2]<<' '<<setw(10)<<I_sd[2][2]<<endl;
	
	// INITIALIZE PARTICLES
    counter=0;
	while(true){
        if(counter==0){
            for(j=0;j<NPAR;j++){
                params[j]=gsl_ran_flat(r,0,30);
            }
        } else { 
            for(j=0;j<NPAR;j++) {
                if(theta[j]<1e-2) params[j]=gsl_ran_flat(r,0,0.01);
                else params[j]=theta[j]*gsl_ran_flat(r,.9,1.1);
            }
        }
		
		solveCTM(params,0,0,out_array_mc);
		solveCTM(params,1,0,out_array_mc_KOA);
		solveCTM(params,0,1,out_array_mc_KOB);
		solveCTM(params,1,1,out_array_mc_KOAB);
		
		model=get_ID(params);
		
		feat_vec[0]=out_array_mc_KOA[0]/out_array_mc[1];
		feat_vec[1]=out_array_mc_KOB[0]/out_array_mc[1];
		feat_vec[2]=out_array_mc_KOAB[0]/out_array_mc[1];
		feat_vec[3]=out_array_mc_KOA[1]/out_array_mc[0];
		feat_vec[4]=out_array_mc_KOB[1]/out_array_mc[0];
		feat_vec[5]=out_array_mc_KOAB[1]/out_array_mc[0];
		feat_vec[6]=out_array_mc_KOA[2]/out_array_mc[2];
		feat_vec[7]=out_array_mc_KOB[2]/out_array_mc[2];
		feat_vec[8]=out_array_mc_KOAB[2]/out_array_mc[2];

		
		likelihood=pow((feat_vec[0]-I[0][0])/I_sd[0][0],2)+
			       pow((feat_vec[1]-I[0][1])/I_sd[0][1],2)+
				   pow((feat_vec[2]-I[0][2])/I_sd[0][2],2)+
				   pow((feat_vec[3]-I[1][0])/I_sd[1][0],2)+
				   pow((feat_vec[4]-I[1][1])/I_sd[1][1],2)+
				   pow((feat_vec[5]-I[1][2])/I_sd[1][2],2)+
				   pow((feat_vec[6]-I[2][0])/I_sd[2][0],2)+
				   pow((feat_vec[7]-I[2][1])/I_sd[2][1],2)+
				   pow((feat_vec[8]-I[2][2])/I_sd[2][2],2);

		// Now since this sample has been calculated from the prior distribution 
		// the weights are just given by the likelihood.

        counter++;

        if(likelihood<likelihood_max){
            for(j=0;j<NPAR;j++){
                theta[j]=params[j];
            }
            likelihood_max=likelihood;
            if(counter>30000){
                if(likelihood<9){
                    for(j=0;j<NPAR;j++) out_file<<params[j]<<' ';
                    out_file<<model<<' '<<likelihood<<' '<<endl;
                }
                likelihood_max=1e10;
                counter=0;
            }
            
        }

        


	}


	return 0;
}
// ---------------------------------------------------------------------

