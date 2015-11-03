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

### Utils functions
source("sherif_test_FCT.R")


######################
###                ###
###   PARAMETERS   ###
###                ###
######################


### Set simulation parameters
###
param.simul <- list(horizon = 90,  # <-- when the simulations stop
					mc_iter = 3,  # <-- Monte Carlo iterations
					calc_WIW_Re = 0,  # <-- if matrix of Who Infects Who is calculated (WARNING: SLOWS DOWN A LOT)
					timeStepTauLeap = 0.2, # <-- size of the tau-leap approximation
					seed = 1234,
					timeIdxGI = 25, # <-- time when GI is retrieved from simulations
					silentMode = 1)  # <-- silent = no debug messages displayed while running

### Set model parameters 
### (see documentation for input parameter definitions)
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
					betaType = "reduced",  # <-- How parameters 'beta_xx' should be interpreted: "standard" or "reduced" 
					beta_IS = 0.2,
					beta_IwS = 1,
					beta_ISw = 1,
					beta_IwSw = 1,
					beta_FS = 0.2,
					beta_FSw = 1,
					beta_HSw = 0.4,
					beta_timedep = "param_beta_timedep.csv",  # <-- file specifying time-dependence of beta parameters
					overwrite_beta_timedep=c("beta_IS_tstart_vec0","beta_IS_newval_vec2"),
					overwrite_beta_timedep_val=c(5,0.222)
					)
	
### Spatial parameters:
###
param.spatial <- loadSpatialParam()
param.model.sp <- replicate.contact.rates(nLocations = param.spatial[["nLocations"]], param.model)
param.model.sp <- c(param.model.sp, loadParamMigration("gravity_cst.csv"))


######################
###                ###
###  SIMULATIONS   ###
###                ###
######################


### Runs the SHERIF model with specified parameters
### (one location, i.e. no spatial structure)
###
x <- rcpp_sherif(paramsSimul = param.simul, 
				 paramsModel = param.model)

### Runs the *spatial* SHERIF model with specified parameters
###
x.sp <- rcpp_sherif_spatial(paramsSimul = param.simul, 
                            paramsModel = param.model.sp,
                            paramsSpatial = param.spatial)



######################
###                ###
###     PLOTS      ###
###                ###
######################

### Merely plots incidence of first Monte-Carlo iteration
### to quickly check if something went *really* wrong...
###
pdf("sherif_test.pdf")
# non-spatial:
plot(x=x$time[[1]], y=x$cumIncidence[[1]],typ="s",main="Incidence (non-spatial model)")

b <- x$check_beta_IS
b2<- x$check_beta_ISw
b3<- x$check_beta_FS
b4<- x$check_beta_FSw
par(mfrow=c(2,2))
plot(x=1:length(b), y=b,typ="l", main="beta_IS",lwd=6,xlab="days",ylab="",las=1); grid()
plot(x=1:length(b2), y=b2,typ="l", main="beta_ISw",lwd=6,xlab="days",ylab="",las=1); grid()
plot(x=1:length(b3), y=b3,typ="l", main="beta_FS",lwd=6,xlab="days",ylab="",las=1); grid()
plot(x=1:length(b4), y=b4,typ="l", main="beta_FSw",lwd=6,xlab="days",ylab="",las=1); grid()
par(mfrow=c(1,1))
# spatial:
nloc <- length(x.sp)
nloc.p <- ceiling(sqrt(nloc))
par(mfrow=c(nloc.p,nloc.p))
for(i in 1:nloc) plot(x=x.sp[[i]]$time[[1]], 
                      y=x.sp[[i]]$cumIncidence[[1]],
                      main = paste("Incidence at location #",i),
                      typ="s")
dev.off()

if(length(x)>0) 
    message("sherif library seems to be working; check sherif_test.pdf for a quick visual diagnostic...")