//
//  spatialSim.h
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-09-11.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#ifndef __SHERIF__spatialSim__
#define __SHERIF__spatialSim__

#include <stdio.h>
#include "simulator.h"
//#include "mc.h"


class spatialSim{
	
	vector<simulator>		_simulator;
	
	unsigned int			_nLocations;
	Matrix					_distLocations;
	vector<unsigned long>	_popLocations;
	
	vector<double>			_time_firstCase;	// time when first incidence case appeared (for each location)
	
	vector<double>			_gravity_cst;

	unsigned long			_cumMovements;		// Cumulative number of movements of individuals between locations
	
	bool					all_in_R_or_D();
	
public:
	
	spatialSim(){};
	
	spatialSim(unsigned int nLocations,
			   vector<unsigned long> popLocations,
			   vector<double> distLocations,
			   vector<double> migrationParams);
	
	void displayInfo();
	
	void initialize_all_simulators(vector<double>	beta_IS,
								   vector<double>	beta_FS,
								   vector<double>	beta_IwS,
								   vector<double>	beta_ISw,
								   vector<double>	beta_FSw,
								   vector<double>	beta_IwSw,
								   vector<double>	beta_HSw,
								   
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
								   unsigned int nE,
								   unsigned int nI,
								   unsigned int nH,
								   unsigned int nF);
	
	void run_tauLeap_spatial(double horizon,
							 double timestepSize,
							 vector<unsigned long> initI,
							 vector<unsigned long> initIw,
							 vector<unsigned long> initSw,
							 bool calc_WIW_Re);
	
	
	// ==========================
	
	simulator	get_simulator(unsigned int i) {return _simulator[i];}
	
	
	// ==========================
	
	void		set_gravity_cst(vector<double> x){_gravity_cst=x;}
	
	
	// ==========================
	
	void		move_indiv(unsigned int loc_i, unsigned int loc_j, vector<unsigned long> ID);
	
	void		migration_gravity(unsigned int loc_i,
								  unsigned int loc_j,
								  unsigned int state_min,
								  unsigned int state_max,
								  double dt);
	
	
	
	vector<double>		get_time_firstCase();
	
	unsigned long		get_cumMovements(){return _cumMovements;}
};


#endif /* defined(__SHERIF__spatialSim__) */