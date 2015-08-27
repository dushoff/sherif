//
//  simulator.cpp
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "simulator.h"


simulator::simulator(double	beta_IS,
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
					 unsigned int nF
)
{
	/// Constructs a simulator object
	
	_beta_IS = beta_IS;
	_beta_FS = beta_FS;
	_beta_IwS = beta_IwS;
	_beta_ISw = beta_ISw;
	_beta_FSw = beta_FSw;
	_beta_IwSw = beta_IwSw;
	_beta_HSw = beta_HSw;

	_sigma = sigma;
	
	_gamma_H = gamma_H;
	_gamma_Hw = gamma_Hw;
	_gamma_F = gamma_F;
	_gamma_R = gamma_R;
	_h_F = h_F;
	_h_R = h_R;
	_f = f;
	
	_delta = delta;
	_deltaH = deltaH;
	_pH = pH;
	_pHw = pHw;
	
	_popSize = popSize;
	_indiv.resize(popSize);
	
	_nE = nE;
	_nI = nI;
	_nH = nH;
	_nF = nF;
	
	set_eventList();
	
	for(unsigned long i=0; i<popSize; i++)	{
		_indiv[i].create();
		_indiv[i].set_ID(i);
	}
}

void simulator::displayPopulation(bool indivDetails){
	
	cout << "Population census:" <<endl;
	
	cout << "S = " << _count_S << endl;
	cout << "E = " << _count_E << endl;
	cout << "I = " << _count_I << endl;
	cout << "Sw = " << _count_Sw << endl;
	cout << "Ew = " << _count_Ew << endl;
	cout << "Iw = " << _count_Iw << endl;
	
	cout << "H = " << _count_H << endl;
	cout << "F = " << _count_F << endl;
	cout << "R = " << _count_R << endl;
	cout << "D = " << _count_D << endl;
	cout<< endl;
	
	if(indivDetails){
		cout << "ID\t" << "state" <<endl;
		coutline(40);
		for(int i=0; i<_popSize; i++){
			cout<< _indiv[i].get_ID()<<"\t"<<_indiv[i].get_state()<<endl;
		}
	}
	
}


void simulator::displayInfo()
{
	coutline(40);
	cout<<" === SIMULATOR INFO ==="<<endl;
	
	cout << "Population size: "<<_indiv.size()<<endl;
	
	cout << "# E compartments:" << _nE <<endl;
	cout << "# I compartments:" << _nI <<endl;
	cout << "# H compartments:" << _nH <<endl;
	cout << "# F compartments:" << _nF <<endl;
	
	cout << "sigma rates:"; displayVector(_sigma);
	cout << "gamma_H rates:"; displayVector(_gamma_H);
	cout << "gamma_R rates:"; displayVector(_gamma_R);
	cout << "gamma_F rates:"; displayVector(_gamma_F);
	
	cout << "Number of events:" << _count_events_sim << endl;
	cout << "Number of time events: "<<_time.size()<<endl;
	cout << "Last date: "<<_time[_time.size()-1]<<endl;
	cout << "Cumul incidence gen pop: "<<_cumIncidence[_cumIncidence.size()-1]<<endl;
	cout << "Cumul incidence HCW: "<<_cumIncidence_hcw[_cumIncidence_hcw.size()-1]<<endl;
	coutline(40);
}



void simulator::set_eventList(){
	_eventList.push_back("infection_S_by_I");
	_eventList.push_back("infection_S_by_Iw");
	_eventList.push_back("infection_S_by_F");
	_eventList.push_back("infection_Sw_by_I");
	_eventList.push_back("infection_Sw_by_Iw");
	_eventList.push_back("infection_Sw_by_F");
	_eventList.push_back("infection_Sw_by_H");
	_eventList.push_back("progress_E");
	_eventList.push_back("progress_Ew");
	_eventList.push_back("progress_I");
	_eventList.push_back("progress_Iw");
	_eventList.push_back("progress_H");
	_eventList.push_back("progress_F");
	_eventList.push_back("infectOnset");
	_eventList.push_back("infectOnset_HCW");
	_eventList.push_back("hospital");
	_eventList.push_back("hospital_HCW");
	_eventList.push_back("funeral_nonHosp");
	_eventList.push_back("funeral_nonHosp_HCW");
	_eventList.push_back("funeral_Hosp");
	_eventList.push_back("recovery_I");
	_eventList.push_back("recovery_Iw");
	_eventList.push_back("recovery_H");
	_eventList.push_back("deathBuried");
}



unsigned long simulator::findIndivIdx(unsigned long ID){
	
	/// Finds the index position of a given individual ID
	/// in the vector of individuals of the simulation '_indiv'
	
	unsigned long res = _popSize;
	for (unsigned long i=0; i<_popSize; i++) {
		if (_indiv[i].get_ID()==ID) {
			res = i;
			break;
		}
	}
	string errmsg = "ID "+to_string(ID)+ " not found!";
	stopif(res==_popSize,errmsg);
	return res;
}


// Contact rates are constant for now...
double		simulator::_beta_IS_fct(double t){return _beta_IS;}
double		simulator::_beta_FS_fct(double t){return _beta_FS;}
double		simulator::_beta_IwS_fct(double t){return _beta_IwS;}
double		simulator::_beta_ISw_fct(double t){return _beta_ISw;}
double		simulator::_beta_FSw_fct(double t){return _beta_FSw;}
double		simulator::_beta_IwSw_fct(double t){return _beta_IwSw;}
double		simulator::_beta_HSw_fct(double t){return _beta_HSw;}


double simulator::eventRate_infection_S_by_I(){
	/// Infection rate on susceptible by infectious, all in general pop
	
	double t = _currentTime;
	return _beta_IS_fct(t)*_count_I*_count_S/_popSize;
}

double simulator::eventRate_infection_S_by_Iw(){
	/// Infection rate on susceptible by infectious HCW
	
	double t = _currentTime;
	return _beta_IwS_fct(t)*_count_Iw*_count_S/_popSize;
}


double simulator::eventRate_infection_S_by_F(){
	/// Infection rate on susceptible by funeral ceremonies
	
	double t = _currentTime;
	return _beta_FS_fct(t)*_count_F*_count_S/_popSize;
}


double simulator::eventRate_infection_Sw_by_I(){
	/// Infection rate on susceptible by infectious, all in general pop
	
	double t = _currentTime;
	return _beta_ISw_fct(t)*_count_I*_count_Sw/_popSize;
}

double simulator::eventRate_infection_Sw_by_Iw(){
	/// Infection rate on susceptible by infectious, all in general pop
	
	double t = _currentTime;
	return _beta_IwSw_fct(t)*_count_Iw*_count_Sw/_popSize;
}


double simulator::eventRate_infection_Sw_by_F(){
	/// Infection rate on susceptible by infectious, all in general pop
	
	double t = _currentTime;
	return _beta_FSw_fct(t)*_count_F*_count_Sw/_popSize;
}


double simulator::eventRate_infection_Sw_by_H(){
	/// Infection rate on susceptible by infectious, all in general pop
	
	double t = _currentTime;
	return _beta_HSw_fct(t)*_count_H*_count_Sw/_popSize;
}



double simulator::eventRate_progress_E(unsigned int i){
	
	/// Returns the rate to progress from the ith latent stage
	/// (stages from i='1' to '_nE-1')
	
	stopif(i>=_nE || i<=0, "index out of bounds");
	unsigned long NEi = census_state(getState_E(i));
	return _sigma[i-1]*NEi;
}


double simulator::eventRate_progress_Ew(unsigned int i){
	
	/// Returns the rate to progress from the ith latent stage for HCW
	/// (stages from i='1' to '_nE-1')
	
	stopif(i>=_nE || i<=0, "index out of bounds");
	return _sigma[i-1]*census_state(getState_Ew(i));
}



double simulator::eventRate_infectOnset(){
	
	/// Rate on infectiousness onset (E->I)
	
	return _sigma[_nE-1]*census_state(getState_E(_nE));
}



double simulator::eventRate_infectOnset_HCW(){
	
	/// Rate on infectiousness onset (Ew->Iw)
	
	return _sigma[_nE-1]*census_state(getState_Ew(_nE));
}



double simulator::eventRate_progress_I(unsigned int i){
	
	/// Returns the rate to progress from the ith infectious stage
	/// (stages from i='1' to '_nI-1')
	
	stopif(i>=_nI || i<=0, "index out of bounds");
	double gamma = _pH*_gamma_H[i-1] + (1-_pH)*(_delta*_gamma_F[i-1]+(1-_delta)*_gamma_R[i-1]);
	return gamma*census_state(getState_I(i));
}


double simulator::eventRate_progress_Iw(unsigned int i){
	
	/// Returns the rate to progress from the ith infectious stage for HCW
	/// (stages from i='1' to '_nI-1')
	
	stopif(i>=_nI || i<=0, "index out of bounds");
	double gamma = _pHw*_gamma_Hw[i-1] + (1-_pHw)*(_delta*_gamma_F[i-1]+(1-_delta)*_gamma_R[i-1]);
	return gamma*census_state(getState_Iw(i));
}

double simulator::eventRate_hospital(){
	/// Hospitalization event rate for general population
	
	return _pH * _gamma_H[_nI-1] * census_state(getState_I(_nI));
}

double simulator::eventRate_hospital_HCW(){
	/// Hospitalization event rate for general population
	
	return _pHw * _gamma_Hw[_nI-1] * census_state(getState_Iw(_nI));
}


double simulator::eventRate_progress_H(unsigned int i){
	/// Progression rate in hospitalization
	
	double h = _deltaH*_h_F[i-1] + (1-_deltaH)*_h_R[i-1];
	return h * census_state(getState_H(i));
}


double simulator::eventRate_funeral_Hosp(){
	/// Death rate when hospitalized
	
	return _deltaH*_h_F[_nH-1]*census_state(getState_H(_nH));
}

double simulator::eventRate_funeral_nonHosp(){
	/// Death rate when hospitalized
	
	return (1-_pH)*_delta*_gamma_F[_nI-1]*census_state(getState_I(_nI));
}


double simulator::eventRate_funeral_nonHosp_HCW(){
	/// Death rate when hospitalized
	
	return (1-_pHw)*_delta*_gamma_F[_nI-1]*census_state(getState_Iw(_nI));
}

double simulator::eventRate_progress_F(unsigned int i){
	/// Progression rate in funerals
	
	return _f[i-1]*census_state(getState_F(i));
}


double simulator::eventRate_recovery_I(){
	/// Recovery rate of infectious non-hospitalized from general population

	return (1-_delta)*(1-_pH)*_gamma_R[_nI-1]*census_state(getState_I(_nI));
}


double simulator::eventRate_recovery_Iw(){
	/// Recovery rate of infectious non-hospitalized from HCW
	
	return (1-_delta)*(1-_pHw)*_gamma_R[_nI-1]*census_state(getState_Iw(_nI));
}

double simulator::eventRate_recovery_H(){
	/// Recovery rate of hospitalized patients
	
	return (1-_deltaH)*_h_R[_nH-1]*census_state(getState_H(_nH));
}

double simulator::eventRate_deathBuried(){
	/// Burial rate of dead after funerals
	
	return _f[_nF-1]*census_state(getState_F(_nF));
}


double simulator::eventRate(string eventType, unsigned int i=0){
	/// Link event rates functions
	
	double rate = -9E9;
	
	if(eventType=="infection_S_by_I")	rate = eventRate_infection_S_by_I();
	if(eventType=="infection_S_by_Iw")	rate = eventRate_infection_S_by_Iw();
	if(eventType=="infection_S_by_F")	rate = eventRate_infection_S_by_F();
	
	if(eventType=="infection_Sw_by_I")	rate = eventRate_infection_Sw_by_I();
	if(eventType=="infection_Sw_by_Iw")	rate = eventRate_infection_Sw_by_Iw();
	if(eventType=="infection_Sw_by_F")	rate = eventRate_infection_Sw_by_F();
	if(eventType=="infection_Sw_by_H")	rate = eventRate_infection_Sw_by_H();
	
	if(eventType=="progress_E")			rate = eventRate_progress_E(i);
	if(eventType=="progress_Ew")		rate = eventRate_progress_Ew(i);
	if(eventType=="progress_I")			rate = eventRate_progress_I(i);
	if(eventType=="progress_Iw")		rate = eventRate_progress_Iw(i);
	if(eventType=="progress_H")			rate = eventRate_progress_H(i);
	if(eventType=="progress_F")			rate = eventRate_progress_F(i);
	
	if(eventType=="infectOnset")		rate = eventRate_infectOnset();
	if(eventType=="infectOnset_HCW")	rate = eventRate_infectOnset_HCW();

	if(eventType=="hospital")			rate = eventRate_hospital();
	if(eventType=="hospital_HCW")		rate = eventRate_hospital_HCW();

	if(eventType=="funeral_nonHosp_HCW")rate = eventRate_funeral_nonHosp_HCW();
	if(eventType=="funeral_nonHosp")	rate = eventRate_funeral_nonHosp();
	if(eventType=="funeral_Hosp")		rate = eventRate_funeral_Hosp();
	
	if(eventType=="recovery_I")			rate = eventRate_recovery_I();
	if(eventType=="recovery_Iw")		rate = eventRate_recovery_Iw();
	if(eventType=="recovery_H")			rate = eventRate_recovery_H();

	if(eventType=="deathBuried")		rate = eventRate_deathBuried();

	return rate;
}



EventNumberContainer simulator::drawNumberEvents_tauLeap(double timestep)
{
	/// Draws the number of events during a tau-leap
	/// for every event
	
	vector<unsigned long> n_event;
	vector<string> event_label;
	vector<unsigned int> event_sublabel;
	bool isProgressEvent;
	
	for(int i=0;i<_eventList.size();i++){
		
		isProgressEvent = is_in_string(_eventList[i], "progress");
		
		if(!isProgressEvent) {
			n_event.push_back(poisson(eventRate(_eventList[i])*timestep));
			event_label.push_back(_eventList[i]);
			event_sublabel.push_back(0);
		}
		if(isProgressEvent){
			if (_eventList[i]=="progress_E" || _eventList[i]=="progress_Ew") {
				for(int k=1;k<_nE;k++){
					n_event.push_back(poisson(eventRate(_eventList[i],k)*timestep));
					event_label.push_back(_eventList[i]);
					event_sublabel.push_back(k);
				}
			}
			if (_eventList[i]=="progress_I" || _eventList[i]=="progress_Iw") {
				for(int k=1;k<_nI;k++){
					n_event.push_back(poisson(eventRate(_eventList[i],k)*timestep));
					event_label.push_back(_eventList[i]);
					event_sublabel.push_back(k);
				}
			}
			if (_eventList[i]=="progress_H") {
				for(int k=1;k<_nH;k++){
					n_event.push_back(poisson(eventRate(_eventList[i],k)*timestep));
					event_label.push_back(_eventList[i]);
					event_sublabel.push_back(k);
				}
			}
			if (_eventList[i]=="progress_F") {
				for(int k=1;k<_nF;k++){
					n_event.push_back(poisson(eventRate(_eventList[i],k)*timestep));
					event_label.push_back(_eventList[i]);
					event_sublabel.push_back(k);
				}
			}
		}
	}
	
	EventNumberContainer ENC(event_label, event_sublabel, n_event);
	return ENC;
}



void simulator::clean_start()
{
	/// Makes sure the simulation starts clean
	
	_prevalence.clear();
	_time.clear();
	
	_cumIncidence.clear();
	_cumIncidence_hcw.clear();
	_incidence.clear();
	_incidence_hcw.clear();
	

	_count_S_vec.clear();
	_count_Sw_vec.clear();
	_count_E_vec.clear();
	_count_Ew_vec.clear();
	_count_I_vec.clear();
	_count_Iw_vec.clear();
	_count_R_vec.clear();
	_count_D_vec.clear();
	_count_F_vec.clear();
	_fatalCum_genPop_vec.clear();
	_fatalCum_hcw_vec.clear();
	
	_WIW_times.clear();
	_WIW.clear();
	_Reff_final.clear();
	
	_GIbck.clear();
	_GIbck_times.clear();

	for(int i=0;i<_popSize; i++) _indiv[i].create();
}



void simulator::initialize(unsigned long initI,
						   unsigned long initIw,
						   unsigned long initSw)
{
	/// Initialize before running the simulation
	
	clean_start();
	
	// Set initial infectious individuals
	
	for(int i=0; i<initI; i++){
		// Warning: initial infectious individuals
		// are put in I[1] (_not_ exposed/latent stage)
		_indiv[i].set_state(getState_I(1));
		_indiv[i].set_timeDiseaseAcquisition(0.0);
	}
	for(int i=0; i<initIw; i++){
		// Warning: initial infectious individuals
		// are put in Iw[1] (_not_ exposed/latent stage)
		_indiv[initI+i].set_state(getState_Iw(1));
		_indiv[initI+i].set_isHCW(true);
		_indiv[initI+i].set_timeDiseaseAcquisition(0.0);
	}
	
	// Set initial susceptible HCW
	for(int i=0; i<initSw; i++){
		// Warning: initial infectious individuals
		// are put in I[1] (_not_ exposed/latent stage)
		_indiv[initI+initIw+i].set_state(getState_Sw());
		_indiv[initI+initIw+i].set_isHCW(true);
		_indiv[initI+initIw+i].set_timeDiseaseAcquisition(0.0);
	}

	
	// intialize counts
	_count_S = _popSize - initI - initIw - initSw;
	_count_E = 0;
	_count_I = initI;
	_count_H = 0;
	_count_F = 0;
	_count_R = 0;
	_count_D = 0;
	
	_fatalCum_genPop = 0;
	_fatalCum_hcw = 0;

	_count_Sw = initSw;
	_count_Ew = 0;
	_count_Iw = initIw;

	
	_eventCount_IS = 0;
	_eventCount_FS = 0;
	_eventCount_IwS = 0;
	_eventCount_ISw = 0;
	_eventCount_FSw = 0;
	_eventCount_IwSw = 0;
	_eventCount_HSw = 0;

	
	// Initialization
	
	double t=0.0;
	_time.push_back(t);
	
	_currCumInc = initI;
	_currCumInc_hcw = initIw;
	_cumIncidence.push_back(initI);
	_cumIncidence_hcw.push_back(initIw);
	_incidence.push_back(initI);
	_incidence_hcw.push_back(initIw);
	
	_prevalence.push_back(initI+initIw);
	
	_count_S_vec.push_back(_count_S);
	_count_Sw_vec.push_back(_count_Sw);
	_count_D_vec.push_back(_count_D);
	_count_R_vec.push_back(_count_R);
	_count_F_vec.push_back(_count_F);
	
	_fatalCum_genPop_vec.push_back(_fatalCum_genPop);
	_fatalCum_hcw_vec.push_back(_fatalCum_hcw);
	
	initialize_table_state_ID();
	
	_count_events_sim = 0;
	_count_events_sim_ignored = 0;
}



void simulator::initialize_table_state_ID(){
	
	/// Iniatializes '_table_state_id'
	/// (called at the start of the simulation only!)
	
	unsigned ns = getTotalNumberOfStates();
	vector< vector<unsigned long> > x(ns);
	
	for(unsigned int s=0; s<ns; s++){
		// Retrieve all individual IDs of state 's'
		x[s] = census_ID(s);
	}
	_table_state_ID = x;
}



void simulator::update_table_state_ID(unsigned long ID,
									  unsigned int prevState,
									  unsigned int nextState){
	
	/// Update '_table_state_ID' for an indiv changing state
	
	// Remove ID from previous state vector:
	_table_state_ID[prevState] = popElementValue(_table_state_ID[prevState], ID);
	// Add ID to next state vector:
	_table_state_ID[nextState].push_back(ID);
}


bool simulator::identify_infector_infectee(string eventType,
										   double timeEvent){
	/// Identify who gets infected by whom
	/// If the event action was triggered
	/// when there is not an associated individual
	/// (as can happen in leap frog algo)
	/// then return 'false' to signal event
	/// could not be actioned
	
	unsigned int state_infectee;
	unsigned int state_infectee_next; // <- once infectee is infected, what is its next state
	unsigned int state_infector_min;
	unsigned int state_infector_max;
	
	bool eventTypeKnown = false;
	
	if(eventType=="infection_S_by_I") {
		state_infectee		= getState_S();
		state_infectee_next	= getState_E(1);
		state_infector_min	= getState_I(1);
		state_infector_max	= getState_I(_nI);
		eventTypeKnown = true;
		_eventCount_IS++;
	}
	if(eventType=="infection_S_by_Iw") {
		state_infectee		= getState_S();
		state_infectee_next	= getState_E(1);
		state_infector_min	= getState_Iw(1);
		state_infector_max	= getState_Iw(_nI);
		eventTypeKnown = true;
		_eventCount_IwS++;
	}
	if(eventType=="infection_S_by_F") {
		state_infectee		= getState_S();
		state_infectee_next	= getState_E(1);
		state_infector_min	= getState_F(1);
		state_infector_max	= getState_F(_nF);
		eventTypeKnown = true;
		_eventCount_FS++;
	}
	
	if(eventType=="infection_Sw_by_I") {
		state_infectee		= getState_Sw();
		state_infectee_next	= getState_Ew(1);
		state_infector_min	= getState_I(1);
		state_infector_max	= getState_I(_nI);
		eventTypeKnown = true;
		_eventCount_ISw++;
	}
	if(eventType=="infection_Sw_by_Iw") {
		state_infectee		= getState_Sw();
		state_infectee_next	= getState_Ew(1);
		state_infector_min	= getState_Iw(1);
		state_infector_max	= getState_Iw(_nI);
		eventTypeKnown = true;
		_eventCount_IwSw++;
	}
	if(eventType=="infection_Sw_by_F") {
		state_infectee		= getState_Sw();
		state_infectee_next	= getState_Ew(1);
		state_infector_min	= getState_F(1);
		state_infector_max	= getState_F(_nF);
		eventTypeKnown = true;
		_eventCount_FSw++;
	}
	if(eventType=="infection_Sw_by_H") {
		state_infectee		= getState_Sw();
		state_infectee_next	= getState_Ew(1);
		state_infector_min	= getState_H(1);
		state_infector_max	= getState_H(_nH);
		eventTypeKnown = true;
		_eventCount_HSw++;
	}

	// If event type not implemented here, stop.
	string errmsg = "event type "+ eventType + "unknown!";
	stopif(!eventTypeKnown,errmsg);
	
	// retrieve all the susceptible IDs
	vector<unsigned long> x = _table_state_ID[state_infectee]; //census_ID(state_infectee);
	// retrieve all the infectious IDs
	vector<unsigned long> x_I = census_table_state_ID(state_infector_min,state_infector_max);//census_ID(state_infector_min,state_infector_max);
	
	if(x.size()==0) _count_events_sim_ignored++; // cout << "Warning: no infectee exists for event type "+eventType;
	if(x_I.size()==0) _count_events_sim_ignored++; //cout << "Warning:no infector exists for event type "+eventType;

	bool success = (x.size()>0 && x_I.size()>0);
	
	if(success){
		// Pick new infectee randomly
		unsigned long ID_infectee = extractElementRandom(x);
		// Pick new infector randomly
		unsigned long ID_infector = extractElementRandom(x_I);
		
		// Now that we know who infected whom, update all relevant quantities.
		// Infectious status and infector ID:
		set_state(ID_infectee, state_infectee_next);
		set_infectorID(ID_infectee, ID_infector);
		// timing of infections from both view points
		set_timeDiseaseAcquisition(ID_infectee, timeEvent);
		set_timeDiseaseTransmit(ID_infector, timeEvent);
		// generation intervals
		double gi = timeEvent - get_timeDiseaseAcquisition(ID_infector);
		set_GIbck(ID_infectee, gi);
		set_GIfwd(ID_infector, gi);
	}
	return (success);
}


bool simulator::identify_progress(string eventType, unsigned int boxcarIdx){
	/// Identify the individuals that will progress through a box-car compartment
	
	unsigned int currState=0;
	bool eventTypeKnown = false;
	
	if (eventType=="progress_E") {
		stopif(boxcarIdx>=_nE, "wrong state for E progression");
		currState = getState_E(boxcarIdx);
		eventTypeKnown = true;
	}
	if (eventType=="progress_Ew") {
		stopif(boxcarIdx>=_nE, "wrong state for Ew progression");
		currState = getState_Ew(boxcarIdx);
		eventTypeKnown = true;
	}
	if (eventType=="progress_I") {
		stopif(boxcarIdx>=_nI, "wrong state for I progression");
		currState = getState_I(boxcarIdx);
		eventTypeKnown = true;
	}
	if (eventType=="progress_Iw") {
		stopif(boxcarIdx>=_nI, "wrong state for Iw progression");
		currState = getState_Iw(boxcarIdx);
		eventTypeKnown = true;
	}
	if (eventType=="progress_H") {
		stopif(boxcarIdx>=_nH, "wrong state for H progression");
		currState = getState_H(boxcarIdx);
		eventTypeKnown = true;
	}
	if (eventType=="progress_F") {
		stopif(boxcarIdx>=_nF, "wrong state for F progression");
		currState = getState_F(boxcarIdx);
		eventTypeKnown = true;
	}
	
	// If event type not implemented here, stop.
	string errmsg = "event type "+ eventType + "unknown!";
	stopif(!eventTypeKnown,errmsg);
	
	// retrieve all candidate IDs
	vector<unsigned long> x = _table_state_ID[currState]; // census_ID(currState);
	
	// IF EXACT ALGO IMPLEMENTED, ACTIVATE LINES BELOW
	// AS AN INTEGRITY CHECK... BUT NEED TO CHANGE BELOW TOO.
	//	errmsg = "no individual exists for event type "+eventType;
	//	stopif(x.size()==0, errmsg);
	
	if(x.size()>0){
		// if x.size()==0 then do nothing.
		// That's justified for the tau-leap algo
		// where the number of event drawn may be
		// larger than the sub-population affected
		// But, this shouldn't be allowed in an exact algo (not implemented yet)
		
		// Pick one randomly
		unsigned long ID = extractElementRandom(x);
		set_state(ID, currState+1);
	}
	
	return (x.size()>0);
}


string simulator::identify_genuineStateChange(string eventType, double timeEvent){
	/// Identify individuals that genuinely change state(e.g., "E->I", "I->F", etc)
	/// not the state change within a box-car compartment.
	/// Exclude infections (dealt by 'identify_infector_infectee').
	///
	/// Three possible outcomes:
	/// - "failed": an indiv that should be affected by the event does not exist (may happen with tau leap, must not happen with exact algo)
	/// - "success_genPop": identification was successfull and the individual affected is from the general population
	/// - "success_hcw": identification was successfull and the individual affected is a HCW
	
	string res = "failed";
	
	unsigned int currState=0;
	unsigned int nextState=0;
	bool eventTypeKnown = false;
	
	if(eventType=="infectOnset") {
		currState = getState_E(_nE);
		nextState = getState_I(1);
		eventTypeKnown = true;
	}
	if(eventType=="infectOnset_HCW") {
		currState = getState_Ew(_nE);
		nextState = getState_Iw(1);
		eventTypeKnown = true;
	}
	if(eventType=="hospital") {
		currState = getState_I(_nI);
		nextState = getState_H(1);
		eventTypeKnown = true;
	}
	if(eventType=="hospital_HCW") {
		currState = getState_Iw(_nI);
		nextState = getState_H(1);
		eventTypeKnown = true;
	}

	if(eventType=="funeral_nonHosp") {
		currState = getState_I(_nI);
		nextState = getState_F(1);
		eventTypeKnown = true;
	}
	if(eventType=="funeral_nonHosp_HCW") {
		currState = getState_Iw(_nI);
		nextState = getState_F(1);
		eventTypeKnown = true;
	}
	if(eventType=="funeral_Hosp") {
		currState = getState_H(_nH);
		nextState = getState_F(1);
		eventTypeKnown = true;
	}
	
	if(eventType=="recovery_I") {
		currState = getState_I(_nI);
		nextState = getState_R();
		eventTypeKnown = true;
	}
	if(eventType=="recovery_Iw") {
		currState = getState_Iw(_nI);
		nextState = getState_R();
		eventTypeKnown = true;
	}
	if(eventType=="recovery_H") {
		currState = getState_H(_nH);
		nextState = getState_R();
		eventTypeKnown = true;
	}
	
	if(eventType=="deathBuried") {
		currState = getState_F(_nF);
		nextState = getState_D();
		eventTypeKnown = true;
	}
	
	// If event type not implemented here, stop.
	string errmsg = "event type "+ eventType + "unknown!";
	stopif(!eventTypeKnown,errmsg);
	
	// retrieve all candidate IDs
	vector<unsigned long> x =  _table_state_ID[currState]; //census_ID(currState);
	// IF EXACT ALGO IMPLEMENTED, ACTIVATE LINES BELOW
	// AS AN INTEGRITY CHECK... BUT NEED TO CHANGE BELOW TOO.
	//	errmsg = "no individual exists for event type "+eventType;
	//	stopif(x.size()==0, errmsg);
	
	if(x.size()==0) _count_events_sim_ignored++; //cout<< "Warning: no individuals associated to an event type "<<eventType<<endl;
	
	if(x.size()>0){
		// if x.size()==0 then do nothing.
		// That's justified for the tau-leap algo
		// where the number of event drawn may be
		// larger than the sub-population affected
		// But, this shouldn't be allowed in an exact algo (not implemented yet)
		
		// Pick one randomly
		unsigned long ID = extractElementRandom(x);
		set_state(ID, nextState);
		// Determine if from general population or HCW:
		res = is_indiv_hcw(ID)?"success_hcw":"success_genPop";
	}
	return res;
}



void simulator::action_event(string eventType, double timeEvent,
							 unsigned int boxcarIdx=0){
	
	/// Performs updates on individuals given an event type
	
	bool event_is_infection_or_progress = false;
	
	_count_events_sim++;
	
	// === INFECTIONS ===
	
	if (eventType=="infection_S_by_I" || eventType=="infection_S_by_Iw" ||
		eventType=="infection_S_by_F") {
		event_is_infection_or_progress = true;
		
		// First, test if individuals associated to the event
		// were found (not necessary true with tau leap algo).
		// If found, then change the counts.
		// If not, do nothing.
		bool success = identify_infector_infectee(eventType, timeEvent);
		if (success && _count_S>0){
			_count_S--;
			_count_E++;
			_currCumInc++;
		}
	}
	
	if (eventType=="infection_Sw_by_I" || eventType=="infection_Sw_by_Iw" ||
		eventType=="infection_Sw_by_F" || eventType=="infection_Sw_by_H") {
		event_is_infection_or_progress = true;
		
		// First, test if individuals associated to the event
		// were found (not necessary true with tau leap algo).
		// If found, then change the counts.
		// If not, do nothing.
		bool success = identify_infector_infectee(eventType, timeEvent);
		if (success && _count_Sw>0){
			_count_Sw--;
			_count_Ew++;
			_currCumInc_hcw++;
		}
	}
	
	// === PROGRESSION through box-car compartments ===
	
	if (eventType=="progress_E" || eventType=="progress_Ew" ||
		eventType=="progress_I" || eventType=="progress_Iw" ||
		eventType=="progress_H" || eventType=="progress_F") {
		identify_progress(eventType, boxcarIdx);
		event_is_infection_or_progress = true;
	}
	
	if(!event_is_infection_or_progress){
		
		// First, test if individuals associated to the event
		// were found (not necessary true with tau leap algo).
		// If found, then change the counts.
		// If not, do nothing.
		string id_res = identify_genuineStateChange(eventType, timeEvent);
		
		if(!(id_res=="failed")){
			// === INFECTIOUSNESS ONSET ===
			
			if (eventType=="infectOnset" && _count_E>0) {
				_count_E--;
				_count_I++;
			}
			if (eventType=="infectOnset_HCW" && _count_Ew>0) {
				_count_Ew--;
				_count_Iw++;
			}
			
			// === HOSPITALIZATION ===
			
			if (eventType=="hospital" && _count_I>0) {
				_count_I--;
				_count_H++;
			}
			if (eventType=="hospital_HCW" && _count_Iw>0) {
				_count_Iw--;
				_count_H++;
			}
			
			// === FUNERALS ===
			
			// Individuals come directly from "I" compartment (necessarily general population)
			if (eventType=="funeral_nonHosp" && _count_I>0) {
				_count_I--;
				_count_F++;
				_fatalCum_genPop++;
			}
			// Individuals come directly from "Iw" compartment (necessarily HCW)
			if (eventType=="funeral_nonHosp_HCW" && _count_Iw>0) {
				_count_Iw--;
				_count_F++;
				_fatalCum_hcw++;
			}
			// Individuals come from "H" compartment (either gen pop or HCW)
			if (eventType=="funeral_Hosp" && _count_H>0) {
				_count_H--;
				_count_F++;
				if(id_res=="success_genPop") {_fatalCum_genPop++;}
				if(id_res=="success_hcw") {_fatalCum_hcw++;}
			}
			
			// === RECOVERY ===
			
			if (eventType=="recovery_I" && _count_I>0) {
				_count_I--;
				_count_R++;
			}
			if (eventType=="recovery_Iw" && _count_Iw>0) {
				_count_Iw--;
				_count_R++;
			}
			if (eventType=="recovery_H" && _count_H>0) {
				_count_H--;
				_count_R++;
			}
			
			// === DEATH & BURIED ===
			
			if (eventType=="deathBuried" && _count_F>0) {
				_count_F--;
				_count_D++;
			}
		}
	}
}



void simulator::run_tauLeap(double horizon,
							double timestepSize,
							unsigned long initI,
							unsigned long initIw,
							unsigned long initSw,
							bool calc_WIW_Re)
{
	/// Run the epidemic simulation with the
	/// tau-leap Poisson approximation of Gillespie algorithm
	
	initialize(initI, initIw, initSw);
	double t = 0.0;

	while ( t<horizon && !all_in_R_or_D() )
	{
		double					time_event = t+timestepSize;
		EventNumberContainer	ENC = drawNumberEvents_tauLeap(timestepSize);
		
		unsigned long nTotalEvents = ENC._event_label.size();
		
		for (int ev=0; ev<nTotalEvents; ev++) {
			
			string eventType = ENC._event_label[ev];
			bool isProgressEvent = is_in_string(eventType, "progress");
		
			// If not a "progress" event type, then action
			// on the event as many times as the number
			// that was drawn from the Poisson in "drawNumberEvents_tauLeap"
			if(!isProgressEvent)
				for(int k=0;k<ENC._n_event[ev];k++)
					action_event(eventType,time_event);
			
			// If a "progress" event type, then must action for
			// each of the box-car compartments:
			if(isProgressEvent){
				for(int k=0;k<ENC._n_event[ev];k++)
					action_event(eventType,time_event,ENC._event_sublabel[ev]); 
			}

			// integrity check (optional, switch off for performance)
			//check_popSize();
			//check_popSize_I();
		}
		
		// update backward GI
		// (not at all event dates b/c of memory cost)
		int tt = round(t);
		if(fabs(t-tt)<0.0001) {
			update_GIbck(t);
		}

		update_incidences();
		_prevalence.push_back(_count_I+_count_Iw+_count_H);  // update prevalence time series
		update_all_count_vec();   // update counts time series
		
		// update times
		_time.push_back(time_event);
		t += timestepSize;
	} // end while
	
	
	// Calculate time series of case effective reproductive number _Reff_final
	// (_WIW must be calculated until horizon)
	if(calc_WIW_Re) {
		calc_WIW(t);
		calc_Reff_final();
		// Integrity check
		if (_WIW.size()>0){
			unsigned long n_wiw = _WIW[_WIW.size()-1].countNonZeroElements();
			unsigned long n_cuminc =_cumIncidence.back();
			unsigned long n_cuminc_hcw =_cumIncidence_hcw.back();
			stopif(n_wiw != n_cuminc+n_cuminc_hcw, "cumulative incidences (_WIW vs _cumIncidence) not consistent!");
		}
	}
	
	if (_count_events_sim_ignored>0){
		cout << "Warning tau-leap algorithm: "<< endl<< _count_events_sim_ignored;
		cout << " events were ignored out of a total of "<<_count_events_sim<<" events (";
		cout << (double)(_count_events_sim_ignored)/(double)(_count_events_sim)*100.0<<"%)"<<endl;
	}
}




void simulator::calc_WIW(double t)
{
	/// Calculate matrix of 'Who Infected Who' at different time points
	
	Matrix M(_popSize,_popSize);
	M.setAllValues(0.0);
	
	unsigned int s	= getState_S();
	unsigned int sw	= getState_Sw();
	
	for(int j=0; j<_popSize; j++){
		// only individuals already infected:
		if (_indiv[j].get_state()!= s &&
			_indiv[j].get_state()!= sw){
			
			unsigned long ii = _indiv[j].get_infectorID();
			M(ii,j) = _indiv[j].get_GIbck();
		}
	}
	_WIW.push_back(M);
	_WIW_times.push_back(t);
}







void simulator::calc_Reff_final()
{
	/// Calculate the realized Effective Reproductive number
	/// It can be calculated only when WIW matrix is updated
	/// and at horizon of the simulation
	
	if (_WIW.size()>0)
	{
	 
		vector<double> n_2nd_cases;
		vector<double> acq_time;
	 
		for (unsigned long i=0; i<_popSize; i++){
		 if (!is_indiv_susceptible(i))
		 {
			 // counts the number of infectees for everyone
			 unsigned long c_i = _WIW[_WIW.size()-1].countNonZeroElements_line(i);
			 n_2nd_cases.push_back((double)c_i);
			 
			 // retrieve disease acquisition time:
			 acq_time.push_back(_indiv[i].get_timeDiseaseAcquisition());
		 }
		}
		_Reff_final.addColVector(acq_time);
		_Reff_final.addColVector(n_2nd_cases);
	}
}



void simulator::update_GIbck(double t){

	/// Update the backward generation interval matrix _GIbck
	/// at time t

	vector<double> gi;
	vector<double> timeAcqInfectee;
	
	// Scan the all infected individuals
	for (unsigned long i=0; i<_popSize; i++){
		if (!is_indiv_susceptible(i)){
			double ta = _indiv[i].get_timeDiseaseAcquisition();
			if(ta>0){
				gi.push_back(_indiv[i].get_GIbck());
				timeAcqInfectee.push_back(ta);
			}
		}
	}
	
	Matrix M(gi);
	M.addColVector(timeAcqInfectee);
	
	_GIbck.push_back(M);
	_GIbck_times.push_back(t);
}




bool simulator::at_least_one_S_and_I_or_F_or_Iw(){
	
	/// Tests if there is at least one susceptible in general pop
	/// AND [one infectious OR Funeral OR infectious HCW]
	/// individual in the whole population
	
	bool res_S = false;
	bool res_2 = false;
	
	unsigned int state_S		= getState_S();
	unsigned int state_Imin		= getState_I(1);
	unsigned int state_Imax		= getState_I(_nI);
	unsigned int state_Iwmin	= getState_Iw(1);
	unsigned int state_Iwmax	= getState_Iw(_nI);
	unsigned int state_Fmin		= getState_F(1);
	unsigned int state_Fmax		= getState_F(_nF);
	
	for(int i=0;i<_popSize; i++){
		unsigned int istatus = _indiv[i].get_state();
		
		if (istatus==state_S) res_S = true;
		
		if ((istatus>= state_Imin && istatus<= state_Imax) ||
			(istatus>= state_Iwmin && istatus<= state_Iwmax) ||
			(istatus>= state_Fmin && istatus<= state_Fmax))
			res_2 = true;
		
		if(res_S && res_2) break;
	}
	return res_S*res_2;
}


bool simulator::all_in_R_or_D(){
	/// Are all individuals in either R or D compartments?
	
	return (_count_D+_count_R == _popSize);
}


void simulator::update_all_count_vec(){
	/// Update time series of subpopulations counts
	
	_count_S_vec.push_back(_count_S);
	_count_E_vec.push_back(_count_E);
	_count_I_vec.push_back(_count_I);
	_count_Sw_vec.push_back(_count_Sw);
	_count_Ew_vec.push_back(_count_Ew);
	_count_Iw_vec.push_back(_count_Iw);
	_count_H_vec.push_back(_count_H);
	_count_F_vec.push_back(_count_F);
	_count_R_vec.push_back(_count_R);
	_count_D_vec.push_back(_count_D);
	
	_fatalCum_genPop_vec.push_back(_fatalCum_genPop);
	_fatalCum_hcw_vec.push_back(_fatalCum_hcw);
}

void simulator::update_incidences(){
	/// Update all incidences
	/// (once _currCumInc is calculated)
	
	// update cumulative incidence
	_cumIncidence.push_back(_currCumInc);
	_cumIncidence_hcw.push_back(_currCumInc_hcw);

	// update period incidence
	unsigned long n = _cumIncidence.size();
	_incidence.push_back(_cumIncidence[n-1]-_cumIncidence[n-2]);
	_incidence_hcw.push_back(_cumIncidence_hcw[n-1]-_cumIncidence_hcw[n-2]);
}



void simulator::check_popSize(){
	/// integrity check for population size
	unsigned long count = _count_S + _count_E + _count_I + _count_Sw + _count_Ew + _count_Iw + _count_H + _count_F + _count_R + _count_D;
	string errmsg = "Population actual size ("+ to_string(_popSize) + ") and counter ("+to_string(count)+") not consistent";
	stopif(count!=_popSize, errmsg);
}

void simulator::check_popSize_I(){
	/// integrity check for infectious general population size
	unsigned int s_min = getState_I(1);
	unsigned int s_max = getState_I(_nI);
	unsigned long count = census_state(s_min, s_max);
	
	string errmsg = "Population actual size ("+ to_string(count) + ") and counter ("+to_string(_count_I)+") for I not consistent";
	stopif(count!=_count_I, errmsg);
}



unsigned long simulator::census_state(unsigned int a, unsigned int b)
{
	/// counts individuals b/w state a and b (both included)
	
	unsigned long cnt = 0;
	
	for(int i=0;i<_popSize; i++){
		unsigned int istatus = _indiv[i].get_state();
		if (istatus>=a && istatus<=b) cnt++;
	}
	return cnt;
}



unsigned long simulator::census_state(unsigned int a)
{
	/// counts individuals of state a
	
	unsigned long cnt = 0;
	
	for(int i=0;i<_popSize; i++){
		if (_indiv[i].get_state()==a) cnt++;
	}
	return cnt;
}



vector<unsigned long> simulator::census_ID(unsigned int a)
{
	/// Retrieve IDs of all individuals of state 'a'
	
	vector<unsigned long> res;
	
	for(int i=0;i<_popSize; i++){
		if (_indiv[i].get_state()==a)
			res.push_back(_indiv[i].get_ID());
	}
	return res;
}

	
vector<unsigned long> simulator::census_ID(unsigned int a, unsigned int b)
{
	/// Retrieve IDs of all individuals between state 'a'and 'b' (a<=b)
	
	stopif(a>b,"states not in right order!");
	
	vector<unsigned long> res;
	
	for(int i=0;i<_popSize; i++){
		if (_indiv[i].get_state()>=a && _indiv[i].get_state()<=b)
			res.push_back(_indiv[i].get_ID());
	}
	return res;
}

vector<unsigned long> simulator::census_table_state_ID(unsigned int a, unsigned int b){
	
	/// Retrieve IDs of all individuals between state 'a'and 'b' (a<=b)
	/// using _table_state_ID
	
	stopif(a>b,"states not in right order!");
	
	vector<unsigned long> res = _table_state_ID[a];
	vector<unsigned long> tmp;
	
	for(unsigned int i=a+1; i<=b; i++){
		tmp = bindVector(res, _table_state_ID[i]);
		res = tmp;
	}
	return res;
}



bool simulator::is_indiv_hcw(unsigned long ID){

	/// Is this individual a HCW?
	
	bool res = false;
	
	unsigned long ii = findIndivIdx(ID);
	res = _indiv[ii].get_isHCW();
	return res;
}


bool simulator::is_indiv_susceptible(unsigned long ID){
	/// Is this individual susceptible?
	
	unsigned int s		= getState_S();
	unsigned int sw		= getState_Sw();
	unsigned long idx	= findIndivIdx(ID);
	
	bool res = true;
	if( _indiv[idx].get_state()!=s &&
	    _indiv[idx].get_state()!=sw)
		res = false;

	return res;
}


void simulator::set_GIbck(unsigned long IDindiv, double gi)
{
	/// Set backward generation interval to value 'gi' for individual ID 'IDindiv'
	
	unsigned long ii = findIndivIdx(IDindiv);
	_indiv[ii].set_GIbck(gi);
}

void simulator::set_GIfwd(unsigned long IDindiv, double gi)
{
	/// Set backward generation interval to value 'gi' for individual ID 'IDindiv'
	
	unsigned long ii = findIndivIdx(IDindiv);
	_indiv[ii].set_GIfwd_incr(gi);
}



double simulator::get_timeDiseaseAcquisition(unsigned long ID)
{
	double res = -9.99;
	
	unsigned long ii = findIndivIdx(ID);
	res = _indiv[ii].get_timeDiseaseAcquisition();
	return res;
}




void simulator::set_state(unsigned long IDindiv, unsigned int s)
{
	/// Set state 's' for a given individual
	
	unsigned long ii = findIndivIdx(IDindiv);
	unsigned int prevState = _indiv[ii].get_state();
	_indiv[ii].set_state(s);
	
	update_table_state_ID(IDindiv, prevState, s);
}



void simulator::set_infectorID(unsigned long IDinfectee, unsigned long IDinfector)
{
	/// Set the infector ID of a given infectee
	
	unsigned long ii = findIndivIdx(IDinfectee);
	_indiv[ii].set_infectorID(IDinfector);
}


void simulator::set_timeDiseaseAcquisition(unsigned long IDinfectee, double t)
{
	unsigned long ii = findIndivIdx(IDinfectee);
	_indiv[ii].set_timeDiseaseAcquisition(t);
}

void simulator::set_timeDiseaseTransmit(unsigned long IDinfector, double t)
{
	unsigned long ii = findIndivIdx(IDinfector);
	_indiv[ii].set_timeDiseaseTransmit(t);
}



unsigned long simulator::census_hcw(){
	
	/// counts all HCW
	
	unsigned long sum = 0;
	for(int i=0; i<_popSize; i++) if(_indiv[i].get_isHCW()) sum++;
	return sum;
}




void simulator::display_eventCounts(){
	
	coutline(40);
	cout << "Infection events counts:" <<endl;
	
	unsigned long total =	_eventCount_IS+_eventCount_FS+
							_eventCount_IwS+_eventCount_ISw+
							_eventCount_FSw+_eventCount_IwSw+
							_eventCount_HSw;
	
	cout << "I -> S:\t"<<_eventCount_IS<<" ("<< round((double)(_eventCount_IS)/total*100) <<"%)"<<endl;
	cout << "F -> S:\t"<<_eventCount_FS<<" ("<< round((double)(_eventCount_FS)/total*100) <<"%)"<<endl;
	cout << "Iw -> S:\t"<<_eventCount_IwS<<" ("<< round((double)(_eventCount_IwS)/total*100) <<"%)"<<endl;
	cout << "I -> Sw:\t"<<_eventCount_ISw<<" ("<< round((double)(_eventCount_ISw)/total*100) <<"%)"<<endl;
	cout << "Iw -> Sw:\t"<<_eventCount_IwSw<<" ("<< round((double)(_eventCount_IwSw)/total*100) <<"%)"<<endl;
	cout << "F -> Sw:\t"<<_eventCount_FSw<<" ("<< round((double)(_eventCount_FSw)/total*100) <<"%)"<<endl;
	cout << "H -> Sw:\t"<<_eventCount_HSw<<" ("<< round((double)(_eventCount_HSw)/total*100) <<"%)"<<endl;
	coutline(40);
}



void simulator::save_GIbck(string filename){
	/// Save the temporal evolution of
	/// backward GI to a file

	ofstream f(filename);
	
	for (int i=0; i<_popSize;  i++){
		double tmp = _indiv[i].get_timeDiseaseAcquisition();
		
		if(tmp>0){
			f << tmp;
			f << ",";
			f << _indiv[i].get_GIbck() << endl;
		}
	}
}




void simulator::save_GIfwd(string filename){
	/// Save the temporal evolution of
	/// forward GI to a file
	
//	ofstream f(filename);
//	
//	double R0 = _beta/_gamma[0]*_nI;
//	// make sure the max length
//	// of GIfwd is sufficient (for all MC iters!)
//	unsigned long gimaxsize = 20*(1+R0);
//	
//	for (int i=0; i<_popSize;  i++)
//		gimaxsize = max(gimaxsize,_indiv[i].get_GIfwd().size());
//	
//	// Now save all GI (even if there's not as much as 'gimaxsize')
//	
//	for (int i=0; i<_popSize;  i++)
//	{
//		double tmp = _indiv[i].get_timeDiseaseAcquisition();
//		unsigned long ntransm = _indiv[i].get_GIfwd().size();
//		
//		if(ntransm>0){
//			f << tmp <<",";
//			
//			for(int k=0;k<gimaxsize;k++)
//			{
//				if(k<ntransm) f<<_indiv[i].get_GIfwd()[k];
//				if(k>=ntransm) f<<"";
//				
//				if(k<gimaxsize-1) f<<",";
//			}
//			f<<endl;
//		}
//	}
}


void simulator::save_sim_outputs(string filename){
	
	ofstream f(filename);
	stopif(_time.size() != _cumIncidence.size(),"simulation inconsistent");
	
	// headers
	f<<"time,cumIncidence,cumIncidence_hcw,count_I,count_Iw,";
	f<<"count_S,count_Sw,count_E,count_Ew,";
	f<<"count_H,count_F,count_R,count_D,prevalence,fatalCum,fatalCum_hcw";
	
	f<< endl;
	
	// data ( SAME ORDER AS HEADERS!!!)
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t];
		f<<",";
		f<<_cumIncidence[t];
		f<<",";
		f<<_cumIncidence_hcw[t];
		f<<",";
		f<<_count_I_vec[t];
		f<<",";
		f<<_count_Iw_vec[t];
		f<<",";
		f<<_count_S_vec[t];
		f<<",";
		f<<_count_Sw_vec[t];
		f<<",";
		f<<_count_E_vec[t];
		f<<",";
		f<<_count_Ew_vec[t];
		f<<",";
		f<<_count_H_vec[t];
		f<<",";
		f<<_count_F_vec[t];
		f<<",";
		f<<_count_R_vec[t];
		f<<",";
		f<<_count_D_vec[t];
		f<<",";
		f<<_prevalence[t];
		f<<",";
		f<<_fatalCum_genPop_vec[t];
		f<<",";
		f<<_fatalCum_hcw_vec[t];
		//f<<",";
		f<<endl;
	}
	
}


void simulator::save_prevalence(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _prevalence.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_prevalence[t]<<endl;
	}
}


void simulator::save_cumIncidence(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _cumIncidence.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_cumIncidence[t]<<endl;
	}
}


void simulator::save_nS(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _count_S_vec.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_count_S_vec[t]<<endl;
	}
}


void simulator::save_nR(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _count_R_vec.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_count_R_vec[t]<<endl;
	}
}




void simulator::save_Reff_final(string filename)
{
	stopif( (_Reff_final.getNbRows() ==0) && (_WIW.size()>0) ,"Reff_final not calculated --> cannot save associated data!");
	
	if(_Reff_final.getNbRows()>0)	_Reff_final.WriteToFileCSV(filename);
	
	if(_Reff_final.getNbRows()==0)
	{
		// Make sure the file is not completely empty...
		vector<double> tmp(2,0.0);
		_Reff_final.addRowVector(tmp);
		_Reff_final.WriteToFileCSV(filename);
	}
}