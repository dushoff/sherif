//
//  ***** THIS FILE IS FOR TESTING ONLY ****
//
//  main.cpp
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include <iostream>

//#include "simulator.h"
#include "individual.h"
#include "mc.h"
#include "spatialSim.h"

int main(int argc, const char * argv[]) {
	
	system("pwd");
	system("date");
	
	// For performance monitoring
	// - do not delete -
	timeval tim;
	gettimeofday(&tim, NULL);
	double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	// ------------------------------------------
	
	

	
	// Read main simulation parameters from file
	
	string fileparam_simul			= "param_simul.csv";
	string fileparam_model_single	= "param_model_singleLoc.csv";
	string fileparam_model			= "param_model.csv";
	
	// Model parameters
	
	unsigned long popSize		= getParameterFromFile("popSize", fileparam_model);
	
	double latent_mean			= getParameterFromFile("meanDur_latent", fileparam_model);
	
	double infectious_mean_H	= getParameterFromFile("meanDur_infectious_H", fileparam_model);
	double infectious_mean_Hw	= getParameterFromFile("meanDur_infectious_Hw", fileparam_model);
	double infectious_mean_F	= getParameterFromFile("meanDur_infectious_F", fileparam_model);
	double infectious_mean_R	= getParameterFromFile("meanDur_infectious_R", fileparam_model);
	
	double hosp_mean_F			= getParameterFromFile("meanDur_hosp_F", fileparam_model);
	double hosp_mean_R			= getParameterFromFile("meanDur_hosp_R", fileparam_model);
	
	double funeral_mean			= getParameterFromFile("meanDur_funeral", fileparam_model);
	
	int nE						= getParameterFromFile("nE", fileparam_model);
	int nI						= getParameterFromFile("nI", fileparam_model);
	int nH						= getParameterFromFile("nH", fileparam_model);
	int nF						= getParameterFromFile("nF", fileparam_model);
	
	string	betaType			= getParameterFromFile_string("betaType", fileparam_model);
	double	beta_IS				= getParameterFromFile("beta_IS", fileparam_model);
	double	beta_FS				= getParameterFromFile("beta_FS", fileparam_model);
	double	beta_IwS			= getParameterFromFile("beta_IwS", fileparam_model);
	double	beta_ISw			= getParameterFromFile("beta_ISw", fileparam_model);
	double	beta_FSw			= getParameterFromFile("beta_FSw", fileparam_model);
	double	beta_IwSw			= getParameterFromFile("beta_IwSw", fileparam_model);
	double	beta_HSw			= getParameterFromFile("beta_HSw", fileparam_model);
	
	double	delta				= getParameterFromFile("delta", fileparam_model);
	double	deltaH				= getParameterFromFile("deltaH", fileparam_model);
	
	double	pH					= getParameterFromFile("pH", fileparam_model);
	double	pHw					= getParameterFromFile("pHw", fileparam_model);
	
	string	fname_beta_timedep	= getParameterFromFile_string("beta_timedep", fileparam_model);
	
	
	// Overwriting time-dependent beta parameters:
	vector<string> overw_prm;
	overw_prm.push_back("beta_IS_tstart_vec0");
	overw_prm.push_back("beta_IS_newval_vec1");
	vector<double> overw_val;
	overw_val.push_back(15);
	overw_val.push_back(0.0777);
	
	// Simulation parameters
	
	double horizon				= getParameterFromFile("horizon", fileparam_simul);
	unsigned long mc_iter		= getParameterFromFile("mc_iter", fileparam_simul);
	double timeStepTauleap		= getParameterFromFile("timeStepTauLeap", fileparam_simul);
	double GIbck_sampleTime		= getParameterFromFile("timeIdxGI", fileparam_simul);
	unsigned long initI			= getParameterFromFile("init_I1", fileparam_model_single);
	unsigned long initIw		= getParameterFromFile("init_Iw1", fileparam_model_single);
	unsigned long initSw		= getParameterFromFile("init_Sw1", fileparam_model_single);
	bool	calc_WIW_Re			= getParameterFromFile("calc_WIW_Re", fileparam_simul);
	bool	silentMode			= getParameterFromFile("silentMode", fileparam_simul);
	
	
	// === Simulations ===
	
	int jobnum = 1;
	
	bool do_singleLocation = true;
	bool do_multiLocation = false;
	
	
	if(do_singleLocation){
		
		// ==== SINGLE PATCH SIMULATION ====
		
		bool singleLocation = true;
		unsigned long firstID = 0;
		
		simulator SIM = initialize_simulation(betaType,
											  beta_IS,
											  beta_FS,
											  beta_IwS,
											  beta_ISw,
											  beta_FSw,
											  beta_IwSw,
											  beta_HSw,
											  
											  fname_beta_timedep,
											  overw_prm,
											  overw_val,
											  
											  latent_mean,
											  infectious_mean_H,
											  infectious_mean_Hw,
											  infectious_mean_F,
											  infectious_mean_R,
											  hosp_mean_F,
											  hosp_mean_R,
											  funeral_mean,
											  
											  delta,
											  deltaH,
											  pH,
											  pHw,
											  popSize,
											  nE, nI, nH, nF,
											  GIbck_sampleTime,
											  singleLocation,
											  firstID);
		
		
		//
		
		vector<double> tmp = SIM.check_values_beta_IS(100);
		displayVector(tmp);
		tmp = SIM.check_values_beta_IS(100);
		displayVector(tmp);
		
		
		// Choose if execution outputs to files
		// or outputs objects (<- used in R wrappin
	
		bool output_to_files = false;
		
		if(output_to_files)
			MC_run_tauLeap(SIM, mc_iter,horizon,
						   timeStepTauleap,
						   initI, initIw,initSw,
						   jobnum,"fout.out",
						   calc_WIW_Re);
		
		if(!output_to_files){
			int seed = 1234;
			vector<simulator> sim_mc = MC_run_tauLeap_sim(SIM,
														  mc_iter,
														  horizon,
														  timeStepTauleap,
														  initI, initIw,initSw,
														  calc_WIW_Re,
														  seed);
		// DEBUG
			cout <<endl <<" GIs at time "<< sim_mc[0].get_GIbck_times(0)<<" :";
			displayVector(sim_mc[0].get_GIbck_gi(0));
		}
	} //end do_singleLocation
	
	
	if(do_multiLocation){
		
		string file_beta_timedep = "beta_IS_timedep.csv";
		
		// ==== SPATIAL SIMULATION ====
		
		unsigned int nLocation = 3;
		vector<unsigned long> popLocations;
		for(int i =0; i<nLocation; i++) popLocations.push_back(500*(i+1));
		
		Matrix distLoc(3,3);
		distLoc(0,1) = 1.0;
		distLoc(0,2) = 2.3;
		distLoc(1,2) = 2;
		distLoc.display();
		vector<double> distLoc_vec = distLoc.melt();
		
		vector<double> vbeta_IS(nLocation,beta_IS);
		vector<double> vbeta_FS(nLocation,beta_FS);
		vector<double> vbeta_IwS(nLocation,beta_IwS);
		vector<double> vbeta_ISw(nLocation,beta_ISw);
		vector<double> vbeta_FSw(nLocation,beta_FSw);
		vector<double> vbeta_IwSw(nLocation,beta_IwSw);
		vector<double> vbeta_HSw(nLocation,beta_HSw);
		
		vector<double> migrationParams;
		migrationParams.push_back(0.001);
		migrationParams.push_back(5.0);
		
		spatialSim spSim(nLocation, popLocations, distLoc_vec, migrationParams);
		
		
		spSim.initialize_all_simulators(betaType,
										vbeta_IS,
										vbeta_FS,
										vbeta_IwS,
										vbeta_ISw,
										vbeta_FSw,
										vbeta_IwSw,
										vbeta_HSw,
										
										fname_beta_timedep,
										overw_prm,
										overw_val,
										
										latent_mean,
										infectious_mean_H,
										infectious_mean_Hw,
										infectious_mean_F,
										infectious_mean_R,
										hosp_mean_F,
										hosp_mean_R,
										funeral_mean,
										
										delta,
										deltaH,
										pH,
										pHw,
										
										nE, nI, nH, nF,
										GIbck_sampleTime);
		
		// Only first location has initial cases
		vector<unsigned long> vinitI(nLocation,0);
		vector<unsigned long> vinitIw(nLocation,0);
		vector<unsigned long> vinitSw(nLocation,initSw);
		vinitI[0] = initI;
		vinitIw[0] = initIw;
		
		int seed = 1234;
		vector<spatialSim> spSim_mc = MC_run_tauLeap_spatial_sim(spSim,
																 mc_iter,
																 horizon,
																 timeStepTauleap,
																 vinitI, vinitIw,vinitSw,
																 calc_WIW_Re,
																 seed, silentMode);
		
		// DEBUG
		
		for(int i=0;i<spSim_mc.size();i++)
		{
			//		cout<<"First time case MC"<<i<<":";
			//		displayVector(spSim_mc[i].get_time_firstCase());
			//		cout<<"doublecheck:"<<endl;
			//		cout << spSim_mc[i].get_simulator(1).get_time_firstCase()<<endl;
			//		cout<<"cumMovements: "<< spSim_mc[i].get_cumMovements()<<endl;
			
			for(int l=0; l<nLocation; l++){
				cout <<endl<<"loc="<< l <<" GIs at time "<< spSim_mc[i].get_simulator(l).get_GIbck_times(0)<<" :";
				displayVector(spSim_mc[i].get_simulator(l).get_GIbck_gi(0));
			}
		}
		
	} // end do_multiLocation

	
	cout<<endl<<"--- E N D ---"<<endl;
	
	
	// --------------------------------------------------------------
	// COMPUTER TIME MONITORING - do not delete!
	
	gettimeofday(&tim, NULL);
	double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	int minutes = (int)((t2-t1)/60.0);
	double sec = (t2-t1)-minutes*60.0;
	cout << endl << " - - - Computational time : ";
	cout << minutes<<" min "<<sec<<" sec" << endl;
	
	// --------------------------------------------------------------
	
	return 0;
}
