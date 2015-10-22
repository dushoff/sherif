###
### This is  a minimal test to be run
### after library installation to
### check everything is ok
###
### DC: Feel free to edit, but this file 
### should remain pretty simple.
###


### Load the sherif library
###
path.sherif.lib <- "./lib"
library(sherif, lib.loc = path.sherif.lib)


### Set simulation parameters
###
param.simul <- list(horizon = 30,  # <-- when the simulations stop
					mc_iter = 3,  # <-- Monte Carlo iterations
					calc_WIW_Re = 0,  # <-- if matrix of Who Infects Who is calculated (WARNING: SLOWS DOWN A LOT)
					timeStepTauLeap = 0.2, # <-- size of the tau-leap approximation
					seed = 1234,
					timeIdxGI = 25, # <-- time when GI is retrieved from simulations
					silentMode = 1)  # <-- silent = no debug messages displayed while running

### Set model parameters (see documentation for definitions)
###
param.model <- list(popSize = 2000,
					init_I1 = 5,  # <-- initial number of infectious individuals
					init_Iw1 = 2, # <-- initial number of infectious HCW
					init_Sw1 = 50, # <-- initial number of susceptible HCW
					meanDur_latent = 7, # <-- mean latent duration 
					meanDur_infectious_H = 3, 
					meanDur_infectious_Hw = 2,
					meanDur_infectious_F = 7,
					meanDur_infectious_R = 9,
					meanDur_hosp_F = 3,
					meanDur_hosp_R = 5,
					meanDur_funeral = 2,
					nE = 3,
					nI = 3,
					nH = 3, 
					nF = 3,
					delta = 0.6,
					deltaH = 0.5,
					pH = 0.3,
					pHw = 0.9,
					betaType = "standard",  # <-- How parameters 'beta_xx' should be interpreted 
					beta_IS = 0.2,
					beta_FS = 0.2,
					beta_IwS = 0.2,
					beta_ISw = 0.2,
					beta_FSw = 0.2,
					beta_IwSw = 0.1,
					beta_HSw = 0.4
					)
	
### Runs the SHERIF model with specified parameters
###
x <- rcpp_sherif(paramsSimul = param.simul, 
				 paramsModel = param.model)

print(x) 
### Merely plots incidence of first Monte-Carlo iteration
### to quickly check if something went *reslly* wrong...
###
pdf("sherif_test.pdf")
plot(x=x$time[[1]], y=x$cumIncidence[[1]],typ="s")
dev.off()

if(length(x)>0) message("sherif library seems to be working; check sherif_test.pdf for a quick visual diagnostic...")