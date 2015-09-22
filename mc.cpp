//
//  mc.cpp
//  Gillepsie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "mc.h"
#include "globalVar.h"


simulator initialize_simulation(double	beta_IS,
								double	beta_FS,
								double	beta_IwS,
								double	beta_ISw,
								double	beta_FSw,
								double	beta_IwSw,
								double	beta_HSw,
								
								double latent_mean,
								double infectious_mean_H,
								double infectious_mean_Hw,
								double infectious_mean_F,
								double infectious_mean_R,
								double hosp_mean_F,
								double hosp_mean_R,
								double funeral_mean,
								
								double	delta,
								double	deltaH,
								double	pH,
								double	pHw,
								unsigned long popSize,
								unsigned int nE,
								unsigned int nI,
								unsigned int nH,
								unsigned int nF,
								unsigned long firstID){
	
	/// Initialize the simulator before being run.
	/// Use 'human readable' inputs
	
	double sigma0	= 1/latent_mean;
	double gamma0_H	= 1/infectious_mean_H;
	double gamma0_Hw	= 1/infectious_mean_Hw;
	double gamma0_F	= 1/infectious_mean_F;
	double gamma0_R	= 1/infectious_mean_R;
	double h0_R		= 1/hosp_mean_R;
	double h0_F		= 1/hosp_mean_F;
	double f0		= 1/funeral_mean;
	
	vector<double> sigma(nE);
	vector<double> gamma_H(nI);
	vector<double> gamma_Hw(nI);
	vector<double> gamma_F(nI);
	vector<double> gamma_R(nI);
	vector<double> h_F(nH);
	vector<double> h_R(nH);
	vector<double> f(nF);
	
	for (int i=0; i<nE; i++) sigma[i]=sigma0 * nE;
	for (int i=0; i<nI; i++) gamma_H[i]	=gamma0_H * nI;
	for (int i=0; i<nI; i++) gamma_Hw[i]=gamma0_Hw * nI;
	for (int i=0; i<nI; i++) gamma_F[i]	=gamma0_F * nI;
	for (int i=0; i<nI; i++) gamma_R[i]	=gamma0_R * nI;
	for (int i=0; i<nH; i++) h_F[i] = h0_F * nH;
	for (int i=0; i<nH; i++) h_R[i] = h0_R * nH;
	for (int i=0; i<nF; i++) f[i]=f0 * nF;
	
	simulator SIM(beta_IS,
				  beta_FS,
				  beta_IwS,
				  beta_ISw,
				  beta_FSw,
				  beta_IwSw,
				  beta_HSw,
				  sigma,
				  gamma_H, gamma_Hw, gamma_F, gamma_R,
				  h_F, h_R,
				  f,
				  delta, deltaH,
				  pH, pHw,
				  popSize,
				  nE, nI, nH, nF,
				  firstID);
	
	return SIM;
}




vector<simulator> MC_run_tauLeap_sim(simulator S, unsigned long iter_mc,
									 double horizon, double timeStep,
									 unsigned long initI,
									 unsigned long initIw,
									 unsigned long initSw,
									 bool calc_WIW_Re,
									 int seed)
{
	vector<simulator> sim_out;
	
	// Forces the seed to be different for each job
	force_seed_reset(seed);
	
	// Monte-Carlo loop
	for(unsigned long i=0; i<iter_mc; i++){
		cout<<endl<<"MC "<<i+1<<"/"<<iter_mc<<" (tau leap "<<timeStep<<")"<<endl;
		S.run_tauLeap(horizon, timeStep, initI, initIw,initSw, calc_WIW_Re);
		sim_out.push_back(S);

	}
	return sim_out;
}





void MC_run_tauLeap(simulator S, unsigned long iter_mc,
					double horizon, double timeStep,
					unsigned long initI,
					unsigned long initIw,
					unsigned long initSw,
					int jobnum, string fnameout,
					bool calc_WIW_Re)
{
	// Forces the seed to be different for each job
	// (each job is executed independently)
	force_seed_reset(jobnum*7);
	
	// Monte-Carlo loop
	for(unsigned long i=0; i<iter_mc; i++)
	{
		cout<<endl<<"MC "<<i+1<<"/"<<iter_mc<<" (tau leap "<<timeStep<<")"<<endl;
		
		S.run_tauLeap(horizon, timeStep, initI, initIw,initSw, calc_WIW_Re);
		save_all_outputs(S, fnameout, jobnum, iter_mc, i);
		S.displayInfo();
		S.display_eventCounts();
	}
}






vector<spatialSim> MC_run_tauLeap_spatial_sim(spatialSim S,
											  unsigned long iter_mc,
											  double horizon,
											  double timeStep,
											  vector<unsigned long> initI,
											  vector<unsigned long> initIw,
											  vector<unsigned long> initSw,
											  bool calc_WIW_Re,
											  int seed,
											  bool silentMode){
	
	
	vector<spatialSim> sim_out;
	
	// Forces the seed to be different for each job
	force_seed_reset(seed);
	
	// Monte-Carlo loop
	for(unsigned long i=0; i<iter_mc; i++){
		if(!silentMode) cout<<endl<<"MC "<<i+1<<"/"<<iter_mc<<" (tau leap "<<timeStep<<")"<<endl;
		S.run_tauLeap_spatial(horizon, timeStep, initI, initIw,initSw, calc_WIW_Re,silentMode);
		sim_out.push_back(S);
		
		// DEBUG
//		cout<< i<<"-DEBUG time first case:";
//		displayVector(S.get_time_firstCase());
//		cout<<i<<"-DEBUG cumMove: "<<S.get_cumMovements()<<endl;
		// ======
		
	}
	return sim_out;
}







void save_all_outputs(simulator S,
					  string param_set, int jobnum,
					  unsigned long iter_mc,
					  unsigned long current_mc){
	
	/// Save all relevant outputs to files
	
	string post			=  to_string((jobnum-1)*iter_mc+current_mc+1) + ".out";
	string fname		= param_set + "_sim_output";
	string f_Reff_final	= param_set + "__Reff_final";
	
	S.save_sim_outputs(_DIR_OUT + fname + post);
	
	// Realized Reff
	string tmp_Reff_final = _DIR_OUT +f_Reff_final + post;
	S.save_Reff_final(tmp_Reff_final);
}




void save_all_outputs_OLD(simulator S,
						  string param_set, int jobnum,
						  unsigned long iter_mc,
						  unsigned long current_mc)
{
	/// Save all relevant outputs to files
	
	string f_prev	= param_set + "__prev";
	string f_cumInc	= param_set + "__cumInc";
	string f_nS		= param_set + "__nS";
	string f_nR		= param_set + "__nR";
	string f_GIbck	= param_set + "__GIbck";
	string f_GIfwd	= param_set + "__GIfwd";
	string f_Reff_final	= param_set + "__Reff_final";
	
	string post =  to_string((jobnum-1)*iter_mc+current_mc+1) + ".out";
	
	// Prevalence
	string tmp_prev = _DIR_OUT +f_prev + post;
	S.save_prevalence(tmp_prev);
	
	// Cumulative Incidence
	string tmp_cumInc = _DIR_OUT +f_cumInc + post;
	S.save_cumIncidence(tmp_cumInc);
	
	// Susceptible
	string tmp_nS = _DIR_OUT +f_nS + post;
	S.save_nS(tmp_nS);
	
	// Recovered
	string tmp_nR = _DIR_OUT +f_nR + post;
	S.save_nR(tmp_nR);
	
	// Generation interval
	string tmp_gibck = _DIR_OUT +f_GIbck + post;
	S.save_GIbck(tmp_gibck);
	
	string tmp_gifwd = _DIR_OUT +f_GIfwd + post;
	S.save_GIfwd(tmp_gifwd);
	
	// Realized Reff
	string tmp_Reff_final = _DIR_OUT +f_Reff_final + post;
	S.save_Reff_final(tmp_Reff_final);
	
}




void test(simulator S)
{
	for (int i=0;i<10;i++) {
		cout<<i<<"->"<<uniform01()<<endl;
		
	}
}

