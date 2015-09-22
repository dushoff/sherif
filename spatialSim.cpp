//
//  spatialSim.cpp
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-09-11.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "spatialSim.h"
#include "mc.h"


spatialSim::spatialSim(unsigned int nLocations,
					   vector<unsigned long> popLocations,
					   vector<double> distLocations,
					   vector<double> migrationParams){
	
	_nLocations = nLocations;
	_popLocations = popLocations;
	
	// Use vector<double> as input instead of directly
	// a Matrix, because that makes the R wrapping easier.
	Matrix mDistLoc(distLocations,nLocations,nLocations);
	_distLocations = mDistLoc;
	
	_time_firstCase.resize(_nLocations);
	
	_gravity_cst = migrationParams;
	_cumMovements = 0;
	
	//DEBUG:
	displayInfo();
	// ===========
}


void spatialSim::initialize_all_simulators(vector<double>	beta_IS,
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
										   unsigned int nF){
	
	/// Set the SAME model parameters for all locations
	/// (hence assumes same behaviour in all locations)
	
	
	// integrity checks
	string errmsg = "beta vector size not equal to number of locations";
	stopif(beta_IS.size()!=_nLocations, errmsg);
	stopif(beta_FS.size()!=_nLocations, errmsg);
	stopif(beta_IwS.size()!=_nLocations,errmsg);
	stopif(beta_ISw.size()!=_nLocations, errmsg);
	stopif(beta_FSw.size()!=_nLocations, errmsg);
	stopif(beta_IwSw.size()!=_nLocations, errmsg);
	stopif(beta_HSw.size()!=_nLocations, errmsg);
	
	_simulator.resize(_nLocations);
	unsigned long firstID = 0;
	
	
	for(int i=0; i<_nLocations; i++){
		_simulator[i] = initialize_simulation(beta_IS[i], beta_FS[i], beta_IwS[i],
											  beta_ISw[i], beta_FSw[i], beta_IwSw[i],
											  beta_HSw[i],
											  latent_mean,
											  infectious_mean_H,
											  infectious_mean_Hw,
											  infectious_mean_F,
											  infectious_mean_R,
											  hosp_mean_F,
											  hosp_mean_R,
											  funeral_mean,
											  delta,deltaH,
											  pH,pHw,
											  _popLocations[i],
											  nE, nI, nH, nF,
											  firstID);
		firstID += _popLocations[i];
		
		// DEBUG
		//_simulator[i].displayPopulation(true);
		// =====
	}
	
}


void spatialSim::displayInfo(){

	cout <<"Number of locations: "<<_nLocations<<endl;
	cout <<"Initial population in each location:"<<endl;
	displayVector(_popLocations);
	
	cout << "Distance matrix b/w "<< _nLocations <<" locations:"<<endl;
	_distLocations.display();
	
	cout << "Migration model parameters:"<<endl;
	displayVector(_gravity_cst);
}


void spatialSim::run_tauLeap_spatial(double horizon,
									 double timestepSize,
									 vector<unsigned long> initI,
									 vector<unsigned long> initIw,
									 vector<unsigned long> initSw,
									 bool calc_WIW_Re,
									 bool silentMode = true){
	
	/// Run the epidemic simulation for ALL locations with the
	/// tau-leap Poisson approximation of Gillespie algorithm
	
	set_silentMode(silentMode);
	
	double t = 0.0;
	
	// Initialize the infectious population number:
	for(int i=0; i<_nLocations; i++) _simulator[i].initialize(initI[i],
															  initIw[i],
															  initSw[i]);
	_cumMovements = 0;
	
	// DEBUG
//	cout<< "DEBUG init prev:"<<endl;
//	for(int i=0; i<_nLocations; i++){
//		cout<<"loc"<<i;
//		displayVector(_simulator[i].get_prevalence());
//		cout<<"_count_E="<<_simulator[i]._count_E<<endl;
//		cout<<"_count_Ew="<<_simulator[i]._count_Ew<<endl;
//		cout<<"_count_I="<<_simulator[i]._count_I<<endl;
//		cout<<"_count_Iw="<<_simulator[i]._count_Iw<<endl;
//		cout<<"_count_H="<<_simulator[i]._count_H<<endl;
//	}

	
	//=====
	
	
	while ( t<horizon )
	{
		double	time_event = t+timestepSize;
		
		for(int i=0; i<_nLocations; i++) {
			
			if(!_simulator[i].all_in_R_or_D()){
				
				EventNumberContainer	ENC = _simulator[i].drawNumberEvents_tauLeap(timestepSize);
				unsigned long nTotalEvents = ENC._event_label.size();
				
				for (int ev=0; ev<nTotalEvents; ev++) {
					
					string eventType = ENC._event_label[ev];
					bool isProgressEvent = is_in_string(eventType, "progress");
					
					// If not a "progress" event type, then action
					// on the event as many times as the number
					// that was drawn from the Poisson in "drawNumberEvents_tauLeap"
					if(!isProgressEvent)
						for(int k=0;k<ENC._n_event[ev];k++)
							_simulator[i].action_event(eventType,time_event,0);
					
					// If a "progress" event type, then must action for
					// each of the box-car compartments:
					if(isProgressEvent){
						for(int k=0;k<ENC._n_event[ev];k++)
							_simulator[i].action_event(eventType,time_event,ENC._event_sublabel[ev]);
					}
					
					// integrity check (optional, switch off for performance)
					//check_popSize();
					//check_popSize_I();
				}
				
				// update backward GI
				// (not at all event dates b/c of memory cost)
				int tt = round(t);
				if(fabs(t-tt)<0.0001) {
					_simulator[i].update_GIbck(t);
				}
				
				_simulator[i].update_incidences();
				_simulator[i].update_prevalence(); // update prevalence time series
				_simulator[i].update_all_count_vec();   // update counts time series
				_simulator[i].update_time_firstCase(time_event);
				// update times
				_simulator[i]._time.push_back(time_event);
			}
		} // for locations
		
		
		// ---------------------------------------------------------------
		// Migrations between locations
		//
		// WARNING: only exposed individuals migrate:
		unsigned int state_min = _simulator[0].getState_E(1);
		unsigned int state_max = _simulator[0].getState_E(_simulator[0]._nE);
		
		for(unsigned int m=0; m<_nLocations; m++)
			for(unsigned int n=m+1; n<_nLocations; n++)
				migration_gravity(m, n, state_min, state_max, timestepSize);
		// ---------------------------------------------------------------
		
		
		t += timestepSize;
	} // end while
	
	// Calculate time series of case effective reproductive number _Reff_final
	// (_WIW must be calculated until horizon)
	for(int i=0; i<_nLocations; i++) {
		if(calc_WIW_Re) {
			_simulator[i].calc_WIW(t);
			_simulator[i].calc_Reff_final();
			// Integrity check
			if (_simulator[i]._WIW.size()>0){
				unsigned long n_wiw = _simulator[i]._WIW[_simulator[i]._WIW.size()-1].countNonZeroElements();
				unsigned long n_cuminc =_simulator[i]._cumIncidence.back();
				unsigned long n_cuminc_hcw =_simulator[i]._cumIncidence_hcw.back();
				stopif(n_wiw != n_cuminc+n_cuminc_hcw, "cumulative incidences (_WIW vs _cumIncidence) not consistent!");
			}
		}
		
		if (_simulator[i]._count_events_sim_ignored>0 && !_silentMode){
			cout << "Warning tau-leap algorithm: "<< endl<< _simulator[i]._count_events_sim_ignored;
			cout << " events were ignored out of a total of "<<_simulator[i]._count_events_sim<<" events (";
			cout << (double)(_simulator[i]._count_events_sim_ignored)/(double)(_simulator[i]._count_events_sim)*100.0<<"%)";
			cout << " [location #"<<i<<"]"<<endl;
		}
	}
}


void spatialSim::move_indiv(unsigned int loc_i,
							unsigned int loc_j,
							vector<unsigned long> ID){
	
	/// Move individuals whose ID is in 'ID' from location #i to location #j
	
	for(int k=0; k<ID.size(); k++){
		unsigned long id_moved = ID[k];
		individual indiv_moved = _simulator[loc_i].get_individual_ID(id_moved);
		_simulator[loc_j].add_indiv(indiv_moved);
		_simulator[loc_i].remove_indiv(id_moved);
		_popLocations[loc_i]--;
		_popLocations[loc_j]++;
	}
	
	_cumMovements += ID.size();
	
	// DEBUG
//	cout << "Moved "<<ID.size()<<" indiv from location #"<<loc_i<<" to location #"<<loc_j<<endl;
//	displayVector(ID);
	// ======
	
}


void spatialSim::migration_gravity(unsigned int loc_i,
								   unsigned int loc_j,
								   unsigned int state_min,
								   unsigned int state_max,
								   double dt){
	
	/// Migration of EXPOSED and INFECTIOUS individuals
	/// between location #i and #j during period of length 'dt'.
	/// Amount of migration diriven by gravity formulation
	/// k1 * Ni * Nj *exp(-k2 * Dij)
	/// Only individuals of state between 'state_min'
	/// and 'state_max' can migrate.
	
	unsigned long Ni = _popLocations[loc_i];
	unsigned long Nj = _popLocations[loc_j];
	
	double k1 = _gravity_cst[0];
	double k2 = _gravity_cst[1];
	
	double Dij = (loc_i<loc_j)?_distLocations(loc_i,loc_j):_distLocations(loc_j,loc_i);
	
	vector<unsigned long> ID_migr_candidate;
	vector<unsigned long> idx_rand;
	vector<unsigned long> ID_migrant;
	
	// === Migrations from i to j ===
	
	// amount of migration
	double m_ij =k1 * Nj * exp(-k2*Dij);
	unsigned long x_ij = binom( min(1.0, m_ij*dt), Ni);
	
	// select individuals that can migrate based on state
	ID_migr_candidate = _simulator[loc_i].census_table_state_ID(state_min, state_max);
	
	if(ID_migr_candidate.size()>0){
		// draw random sample from this group:
		idx_rand = uniformIntVectorUnique(x_ij, 0, ID_migr_candidate.size()-1);
		// IDs of the individuals selected to migrate:
		for (int k=0; k<idx_rand.size(); k++) ID_migrant.push_back(ID_migr_candidate[idx_rand[k]]);
		
		move_indiv(loc_i, loc_j, ID_migrant);
	}
	
	// === Migrations from j to i ===
	
	ID_migr_candidate.clear();
	idx_rand.clear();
	ID_migrant.clear();
	
	// amount of migration
	m_ij = k1 * Ni * exp(-k2*Dij);
	x_ij = binom(min(1.0, m_ij*dt), Nj);
	
	// select individuals that can migrate based on state
	ID_migr_candidate = _simulator[loc_j].census_table_state_ID(state_min, state_max);
	
	if(ID_migr_candidate.size()>0){
		// draw random sample from this group:
		idx_rand = uniformIntVectorUnique(x_ij, 0, ID_migr_candidate.size()-1);
		// IDs of the individuals selected to migrate:
		for (int k=0; k<idx_rand.size(); k++) ID_migrant.push_back(ID_migr_candidate[idx_rand[k]]);
		
		move_indiv(loc_j, loc_i, ID_migrant);
	}
}


vector<double> spatialSim::get_time_firstCase(){
	/// Retrieve the time when the first case
	/// appeared for all locations
	///
	/// Warning: result may be different
	/// at different simulation times.
	/// (What is likely desired is to retrieve the times at HORIZON)
	
	vector<double> res;
	for(int loc=0; loc<_nLocations; loc++) res.push_back(_simulator[loc].get_time_firstCase());
	return res;
}












