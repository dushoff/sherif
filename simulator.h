//
//  simulator.h
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#ifndef __SHERIF__simulator__
#define __SHERIF__simulator__

#include <stdio.h>
#include "individual.h"
#include "EventNumberContainer.h"
#include "globalVar.h"

class simulator
{
	vector<individual>	_indiv;

	
	// This table stores, at the current simulation time,
	// the IDs for each state.
	// It speeds up the search of the ID of an
	// individual for a given state.
	vector< vector<unsigned long> > _table_state_ID;
	
	unsigned long		_popSize;

	// keeps track of number of individuals in
	// any compartment. If "box-car" compartment
	// count individuals in _any_ of them.
	unsigned long		_count_S;
	unsigned long		_count_Sw;
	unsigned long		_count_E;
	unsigned long		_count_Ew;
	unsigned long		_count_I;
	unsigned long		_count_Iw;
	unsigned long		_count_H;
	unsigned long		_count_F;
	unsigned long		_count_R;
	unsigned long		_count_D;
	
	// Compartment "F" or "D" aggregate deaths from
	// both general population and HCW.
	// These two variables distinguish them:
	unsigned long		_fatalCum_genPop;
	unsigned long		_fatalCum_hcw;
	
	// number of compartments for E, I, H and F
	unsigned int _nE;
	unsigned int _nI;
	unsigned int _nH;
	unsigned int _nF;

	// Contact rates
	double		_beta_IS;
	double		_beta_FS;
	double		_beta_IwS;
	double		_beta_ISw;
	double		_beta_FSw;
	double		_beta_IwSw;
	double		_beta_HSw;
	
	double		_beta_IS_fct(double t);
	double		_beta_FS_fct(double t);
	double		_beta_IwS_fct(double t);
	double		_beta_ISw_fct(double t);
	double		_beta_FSw_fct(double t);
	double		_beta_IwSw_fct(double t);
	double		_beta_HSw_fct(double t);
	
	
	vector<double>		_sigma;		// inv. mean latent duration
	vector<double>		_gamma_H;	// inv. mean duration infectious when hospitalization follows
	vector<double>		_gamma_Hw;	// inv. mean duration infectious when hospitalization follows for HCW
	vector<double>		_gamma_F;	// inv. mean duration infectious when death follows
	vector<double>		_gamma_R;	// inv. mean duration infectious when recovery follows
	vector<double>		_h_F;		// inv. mean hospitalization duration when death follows
	vector<double>		_h_R;		// inv. mean hospitalization duration when recovery follows
	vector<double>		_f;			// inv. mean funeral duration
	
	double		_pH;	// proportion of infectious who get hospitalized
	double		_pHw;	// proportion of infectious HCW who get hospitalized
	double		_delta;	// proportion infectious who die when NOT hospitalized
	double		_deltaH;// proportion infectious who die when hospitalized
	
	
	// events times
	double					_currentTime; // updated in the main loop in 'run' functions
	vector<double>			_time;
	unsigned long			_count_events_sim; // counts the total number of events during one simulation
	unsigned long			_count_events_sim_ignored; // counts the total number of events that were ignored (in tau-leap algo)
	
	// time series
	vector<unsigned long>	_prevalence;
	unsigned long			_currCumInc;
	unsigned long			_currCumInc_hcw;
	vector<unsigned long>	_cumIncidence;
	vector<unsigned long>	_cumIncidence_hcw;
	vector<unsigned long>	_incidence;
	vector<unsigned long>	_incidence_hcw;
	
	vector<unsigned long>	_count_S_vec;	// keep track of the number of susceptible
	vector<unsigned long>	_count_Sw_vec;	// keep track of the number of susceptible HCW
	vector<unsigned long>	_count_E_vec;	// keep track of the number of exposed
	vector<unsigned long>	_count_Ew_vec;	// keep track of the number of exposed HCW
	vector<unsigned long>	_count_I_vec;	// keep track of the number of infectious
	vector<unsigned long>	_count_Iw_vec;	// keep track of the number of infectious HCW
	vector<unsigned long>	_count_H_vec;	// keep track of the number of hospitalized
	vector<unsigned long>	_count_R_vec;	// keep track of the number of recovered
	vector<unsigned long>	_count_F_vec;	// keep track of the number of dead not buried
	vector<unsigned long>	_count_D_vec;	// keep track of the number of dead and buried
	
	vector<unsigned long>	_fatalCum_genPop_vec;
	vector<unsigned long>	_fatalCum_hcw_vec;
	


	unsigned long		findIndivIdx(unsigned long ID);
	
	// Status indices
	// S -> 0
	// Sw -> 1
	// E[k] -> k+1
	// etc...
	
	unsigned int getState_S() {return 0;}
	unsigned int getState_Sw() {return 1;}
	unsigned int getState_E(unsigned int k) {stopif(k>_nE,""); return k+1;}
	unsigned int getState_Ew(unsigned int k) {stopif(k>_nE,"");return _nE+k+1;}
	unsigned int getState_I(unsigned int k) {stopif(k>_nI,"");return 2*_nE+k+1;}
	unsigned int getState_Iw(unsigned int k) {stopif(k>_nI,"");return 2*_nE+_nI+k+1;}
	unsigned int getState_H(unsigned int k) {stopif(k>_nH,"");return 2*(_nE+_nI)+k+1;}
	unsigned int getState_F(unsigned int k) {stopif(k>_nF,"");return 2*(_nE+_nI)+_nH+k+1;}
	unsigned int getState_R() {return 2*(_nE+_nI)+_nH+_nF+2;}
	unsigned int getState_D() {return 2*(_nE+_nI)+_nH+_nF+3;}
	
	unsigned int getTotalNumberOfStates(){return 2*(_nE+_nI)+_nH+_nF+4;}
	
	// Realized effective reproductive number at horizon.
	// 1st column: Disease acquisition date for infector
	// 2nd column: number of secondary cases generated by that infector
	// Number of row = population size
	Matrix				_Reff_final;
		
	// Matrix of 'Who Infected Who' at different time points
	// element WIW(i,j)=0 if indiv#i has NOT infected indiv#j
	// element WIW(i,j)=g if indiv#i infected indiv#j with generation interval=g
	vector<Matrix>		_WIW;
	
	// times when _WIW was recorded
	// (to avoid calculating it at each step,
	//  which would be memory expensive)
	vector<double>		_WIW_times;
	
	// "Backward" Generation Intervals
	// First column : generation interval
	// Second column : time of disease acquisition of the infectee
	// There is a vector of GIbck because it can
	// be calculated at several calendar times
	vector<Matrix>		_GIbck;
	// Calendar times when _GIbck is calculated
	vector<double>		_GIbck_times;
	
	
	void			set_GIbck(unsigned long IDindiv, double gi);
	void			set_GIfwd(unsigned long IDindiv, double gi);
	
	void			set_state(unsigned long IDindiv, unsigned int s);
	void			set_infectorID(unsigned long IDinfectee, unsigned long IDinfector);
	
	void			set_timeDiseaseAcquisition(unsigned long IDinfectee, double t);
	void			set_timeDiseaseTransmit(unsigned long IDinfector, double t);
	
	void			calc_WIW(double t);
	void			calc_Reff_final();

	void			update_GIbck(double t);
	
	// Simulation functions:

	void			initialize(unsigned long initSw,
							   unsigned long initIw,
							   unsigned long initI);
	
	void			initialize_table_state_ID();
	void			update_table_state_ID(unsigned long ID,
										  unsigned int prevState,
										  unsigned int nextState);
	
	
	// List of all possible events.
	vector<string>	_eventList;
	void			set_eventList();

	// Infection Events counters
	unsigned long		_eventCount_IS;
	unsigned long		_eventCount_FS;
	unsigned long		_eventCount_IwS;
	unsigned long		_eventCount_ISw;
	unsigned long		_eventCount_FSw;
	unsigned long		_eventCount_IwSw;
	unsigned long		_eventCount_HSw;
	
	// Event rates
	
	double			eventRate(string eventType,unsigned int i);
	
	double			eventRate_infection_S_by_I();
	double			eventRate_infection_S_by_Iw();
	double			eventRate_infection_S_by_F();
	
	double			eventRate_infection_Sw_by_I();
	double			eventRate_infection_Sw_by_Iw();
	double			eventRate_infection_Sw_by_F();
	double			eventRate_infection_Sw_by_H();
	
	double			eventRate_progress_E(unsigned int i);
	double			eventRate_progress_Ew(unsigned int i);
	
	double			eventRate_infectOnset();
	double			eventRate_infectOnset_HCW();
	
	double			eventRate_progress_I(unsigned int i);
	double			eventRate_progress_Iw(unsigned int i);
	
	double			eventRate_hospital();
	double			eventRate_hospital_HCW();
	
	double			eventRate_progress_H(unsigned int i);
	
	double			eventRate_funeral_nonHosp();
	double			eventRate_funeral_Hosp();
	double			eventRate_funeral_nonHosp_HCW();
	
	double			eventRate_progress_F(unsigned int i);
	
	double			eventRate_recovery_I();
	double			eventRate_recovery_Iw();
	double			eventRate_recovery_H();
	
	double			eventRate_deathBuried();
	
	EventNumberContainer	drawNumberEvents_tauLeap(double timestep);
	
	bool			identify_infector_infectee(string eventType, double timeEvent);
	bool			identify_progress(string eventType, unsigned int boxcarIdx);
	string			identify_genuineStateChange(string eventType, double timeEvent);
	// Note: 'identify_genuineStateChange' is of string type, not bool, because
	// there are more than 2 outcomes only:
	// - failure to identify
	// - identification successful and affects individual from general population
	// - identification successful and affects individual from HCW
	
	void			action_event(string eventType, double timeEvent, unsigned int boxcarIdx);
	
	void			update_all_count_vec();
	void			update_incidences();
	
	void			check_popSize();
	void			check_popSize_I();
	

	
public:
	
	simulator(){};
	
	simulator(double	beta_IS,
			  double	beta_FS,
			  double	beta_IwS,
			  double	beta_ISw,
			  double	beta_FSw,
			  double	beta_IwSw,
			  double	beta_HSw,
			  vector<double>sigma,
			  vector<double>gamma_H,
			  vector<double>gamma_Hw,
			  vector<double>gamma_F,
			  vector<double>gamma_R,
			  vector<double>h_F,
			  vector<double>h_R,
			  vector<double>f,
			  double delta,
			  double deltaH,
			  double pH,
			  double pHw,
			  unsigned long popSize,
			  unsigned int nE,
			  unsigned int nI,
			  unsigned int nH,
			  unsigned int nF);
	
	
	
	
	// ===== GET FUNCTIONS =====
	
	unsigned long			get_popSize(){return _popSize;}
	
	vector<unsigned long>	get_cumIncidence(){return _cumIncidence;}
	vector<unsigned long>	get_cumIncidence_hcw(){return _cumIncidence_hcw;}
	vector<unsigned long>	get_incidence(){return _incidence;}
	vector<unsigned long>	get_incidence_hcw(){return _incidence_hcw;}
	
	vector<unsigned long>	get_prevalence(){return _prevalence;}
	
	vector<unsigned long>	get_fatal_incid_all(){return _count_F_vec;} // <-- Incidence of fatalities for all
	vector<unsigned long>	get_fatal_cum_genPop(){return _fatalCum_genPop_vec;}  // <-- cumulative count of fatalities in general population
	vector<unsigned long>	get_fatal_cum_hcw(){return _fatalCum_hcw_vec;}  // <-- cumulative count of fatalities among HCW
	
	vector<double>			get_time(){return _time;}
	
	vector<Matrix>			get_GIbck() {return _GIbck;}
	vector<double>			get_GIbck_times() {return _GIbck_times;}
	Matrix					get_GIbck(unsigned int i) {return _GIbck[i];}
	double					get_GIbck_times(unsigned int i) {return _GIbck_times[i];}
	vector<double>			get_GIbck_gi(unsigned int i) {return _GIbck[i].extractColumn(0);}
	
	
	// ===== Simulation FUNCTIONS =====
	
	void clean_start();
	
	void run(double horizon,
			 unsigned long initI,
			 unsigned long initIw,
			 unsigned long initSw,
			 bool calc_WIW_Re,
			 bool doExact,
			 double timeStep);
	
	void run_tauLeap(double horizon,
					 double timestepSize,
					 unsigned long initI,
					 unsigned long initIw,
					 unsigned long initSw,
					 bool calc_WIW_Re);
	
	
	
	// ===== OTHER FUNCTIONS =====
	
	unsigned long	census_state(unsigned int a, unsigned int b); // counts individuals b/w infectious status a and b
	unsigned long	census_state(unsigned int a); // counts individuals of infectious status a
	unsigned long	census_hcw();
	
	bool	at_least_one_S_and_I_or_F_or_Iw();
	bool	all_in_R_or_D();

	
	vector<unsigned long> census_ID(unsigned int a); // Retrieve IDs of all individuals of state 'a'
	vector<unsigned long> census_ID(unsigned int a, unsigned int b); // Retrieve IDs of all individuals of state b/w 'a' and 'b'
	vector<unsigned long> census_table_state_ID(unsigned int a, unsigned int b); // Retrieve IDs of all individuals of state b/w 'a' and 'b' using _table_state_ID
	
	bool	is_indiv_hcw(unsigned long ID);
	bool	is_indiv_susceptible(unsigned long ID);
	
	double	get_timeDiseaseAcquisition(unsigned long ID);
	
	
	// Save time series to file
	
	void save_sim_outputs(string filename);
	
	void save_prevalence(string filename);
	void save_cumIncidence(string filename);
	void save_nS(string filename);
	void save_nR(string filename);

	void save_Reff_final(string filename);
	void save_GIbck(string filename);
	void save_GIfwd(string filename);
	
	void displayInfo();
	void displayPopulation(bool indivDetails);
	void display_eventCounts();
	
	void test(){displayVector(_count_F_vec);}
};



#endif /* defined(__SHERIF__simulator__) */
