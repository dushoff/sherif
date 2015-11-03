//
//  mc.h
//  Gillepsie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#ifndef __Gillepsie_SEmInR__mc__
#define __Gillepsie_SEmInR__mc__

#include <stdio.h>
#include "simulator.h"
#include "spatialSim.h"

simulator initialize_simulation(string	betaType,
								double	beta_IS,
								double	beta_FS,
								double	beta_IwS,
								double	beta_ISw,
								double	beta_FSw,
								double	beta_IwSw,
								double	beta_HSw,
								
								string	filename_beta_timedep,
								vector<string> overwrite_beta_timedep,
								vector<double> overwrite_value,
								
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
								double GIbck_sampleTime,
								bool singleLocation,
								unsigned long firstID=0
								);


void MC_run_tauLeap(simulator S,
					unsigned long iter_mc,
					double horizon,
					double timeStep,
					unsigned long initI,
					unsigned long initIw,
					unsigned long initSw,
					int jobnum,
					string fnameout,
					bool calc_WIW_Re);


vector<simulator> MC_run_tauLeap_sim(simulator S,
									 unsigned long iter_mc,
									 double horizon,
									 double timeStep,
									 unsigned long initI,
									 unsigned long initIw,
									 unsigned long initSw,
									 bool calc_WIW_Re,
									 int seed);


vector<spatialSim> MC_run_tauLeap_spatial_sim(spatialSim S,
											  unsigned long iter_mc,
											  double horizon,
											  double timeStep,
											  vector<unsigned long> initI,
											  vector<unsigned long> initIw,
											  vector<unsigned long> initSw,
											  bool calc_WIW_Re,
											  int seed,
											  bool silentMode);




void save_all_outputs(simulator S,
					  string param_set, int jobnum,
					  unsigned long iter_mc,
					  unsigned long current_mc);



void save_all_outputs_OLD(simulator S,
						  string param_set, int jobnum,
						  unsigned long iter_mc,
						  unsigned long current_mc);

void test(simulator S);

#endif /* defined(__Gillepsie_SEmInR__mc__) */
