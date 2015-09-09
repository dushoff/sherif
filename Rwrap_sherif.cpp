/// =================================================================
///
///   WRAPPING FOR R
///
///   Created 2015-08-02 by David Champredon
///
/// =================================================================

#include "simulator.h"
#include "individual.h"
#include "mc.h"

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

List rcpp_sherif(List paramsSimul, List paramsModel) {
	
	// Unpack all model parameters
	unsigned long popSize		= paramsModel["popSize"];
	unsigned long initI			= paramsModel["init_I1"];
	unsigned long initIw		= paramsModel["init_Iw1"];
	unsigned long initSw		= paramsModel["init_Sw1"];
	double latent_mean			= paramsModel["meanDur_latent"];
	double infectious_mean_H	= paramsModel["meanDur_infectious_H"];
	double infectious_mean_Hw	= paramsModel["meanDur_infectious_Hw"];
	double infectious_mean_F	= paramsModel["meanDur_infectious_F"];
	double infectious_mean_R	= paramsModel["meanDur_infectious_R"];
	double hosp_mean_F			= paramsModel["meanDur_hosp_F"];
	double hosp_mean_R			= paramsModel["meanDur_hosp_R"];
	double funeral_mean			= paramsModel["meanDur_funeral"];
	int nE						= paramsModel["nE"];
	int nI						= paramsModel["nI"];
	int nH						= paramsModel["nH"];
	int nF						= paramsModel["nF"];
	double	beta_IS				= paramsModel["beta_IS"];
	double	beta_FS				= paramsModel["beta_FS"];
	double	beta_IwS			= paramsModel["beta_IwS"];
	double	beta_ISw			= paramsModel["beta_ISw"];
	double	beta_FSw			= paramsModel["beta_FSw"];
	double	beta_IwSw			= paramsModel["beta_IwSw"];
	double	beta_HSw			= paramsModel["beta_HSw"];
	double	delta				= paramsModel["delta"];
	double	deltaH				= paramsModel["deltaH"];
	double	pH					= paramsModel["pH"];
	double	pHw					= paramsModel["pHw"];
	
	// Unpack all simulation parameters
	unsigned long mc_iter		= paramsSimul["mc_iter"];
	double	horizon				= paramsSimul["horizon"];
	double	timeStepTauleap		= paramsSimul["timeStepTauLeap"];
	bool	calc_WIW_Re			= paramsSimul["calc_WIW_Re"];
	int		seed				= paramsSimul["seed"];
	unsigned int timeIdxGI		= paramsSimul["timeIdxGI"];
	
	// === Simulations ===
	
	simulator SIM = initialize_simulation(beta_IS,
										  beta_FS,
										  beta_IwS,
										  beta_ISw,
										  beta_FSw,
										  beta_IwSw,
										  beta_HSw,
										  
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
										  nE, nI, nH, nF);
	
	vector<simulator> sim_mc = MC_run_tauLeap_sim(SIM,
												  mc_iter,
												  horizon,
												  timeStepTauleap,
												  initI,
												  initIw,
												  initSw,
												  calc_WIW_Re,
												  seed);
	

	List timeSim(mc_iter);
	List cumInc(mc_iter);
	List cumInc_hcw(mc_iter);
	List inc(mc_iter);
	List inc_hcw(mc_iter);
	List deaths(mc_iter);
	List deaths_hcw(mc_iter);
	List deaths_all(mc_iter);
	List GIbck_time(mc_iter);
	List GIbck_gi(mc_iter);
	List Reff_timeAcq(mc_iter);
	List Reff_n2ndCases(mc_iter);
	
	
	for (int i=0; i<mc_iter; i++) {
		// times of events
		timeSim[i]		= sim_mc[i].get_time();

		// Cumulative incidence:
		cumInc[i]		= sim_mc[i].get_cumIncidence();
		cumInc_hcw[i]	= sim_mc[i].get_cumIncidence_hcw();
		
		inc[i]			= sim_mc[i].get_incidence();
		inc_hcw[i]		= sim_mc[i].get_incidence_hcw();
		
		// Deaths:
		deaths[i]		= sim_mc[i].get_fatal_cum_genPop();
		deaths_hcw[i]	= sim_mc[i].get_fatal_cum_hcw();
		deaths_all[i]	= sim_mc[i].get_fatal_incid_all();
		
		// Times when Backward generation interval is requested
		// and vector of all GI at that time
		GIbck_time[i]	= sim_mc[i].get_GIbck_times(timeIdxGI);
		GIbck_gi[i]		= sim_mc[i].get_GIbck_gi(timeIdxGI);
		
		// Effective reproductive number
		// (realized number of secondary cases for every individual)
		Reff_timeAcq[i]		= sim_mc[i].get_Reff_final_timeAcq();
		Reff_n2ndCases[i]	= sim_mc[i].get_Reff_final_n2ndCases();
		
	}
	
	return List::create(Named("time") = timeSim,
						Named("cumIncidence") = cumInc,
						Named("cumIncidence_hcw") = cumInc_hcw,
						Named("incidence") = inc,
						Named("incidence_hcw") = inc_hcw,
						Named("deathsCum") = deaths,
						Named("deathsCum_hcw") = deaths_hcw,
						Named("deathsIncid_all") = deaths_all,
						Named("GIbck_time") = GIbck_time,
						Named("GIbck_gi") = GIbck_gi,
						Named("Reff_timeAcq") = Reff_timeAcq,
						Named("Reff_n2ndCases") = Reff_n2ndCases
						);
	
}