####################################################
####################################################
####    
###    RUNS SPATIAL SHERIF MODEL WITHIN R 
###
###    Created 2015-00-13 by David Champredon
###
####################################################
####################################################

source("format_NIH.R")

# args <- commandArgs(trailingOnly = TRUE)

t0 <- as.numeric(Sys.time())

### R package for SHERIF is installed locally:
###
path.sherif.lib <- "./lib"

library(sherif, lib.loc = path.sherif.lib)

### Detect the number of CPU cores
###
library(parallel) 
ncpus <- detectCores(all.tests = TRUE, logical = FALSE)

source("run_sherif_spatial_FCT.R")
source("./Calibration/loadParam_FCT.R")

### Load simulation parameters
###
prm.simul <- loadParam("param_simul.csv")
prm.model <- loadParam("param_model.csv")


### Load spatial parameters
###
prm.spatial <- loadSpatialParam()
prm.model <- replicate.contact.rates(nLocations = prm.spatial[["nLocations"]], prm.model)
prm.model=c(prm.model,loadParamMigration("gravity_cst.csv"))

### Run the SHERIF simulations
###
sim <- run.sherif.spatial.parallel(prm.simul = prm.simul, 
                                   prm.model = prm.model,
                                   prm.spatial = prm.spatial,
                                   ncpus = ncpus, 
                                   path.sherif.lib = path.sherif.lib)




### Save output in NIH's format
proba.inc <- seq(0,0.12, by=0.001)
proba.cuminc <- seq(0,1, by=0.01)
prediction.date <- read.table("prediction_date.csv")[1,1]
bucket.size <- 7


for(loc in 1:prm.spatial[["nLocations"]]){

	filesuffix <- paste0("REGION_",loc)
	
	NIH.format.sherif.output(sim[[loc]], 
							 popsize = prm.spatial[["popLocations"]][[loc]],
							 proba.inc,
							 proba.cuminc,
							 prediction.date,
							 bucket.size,
							 filesuffix)
}

### PLOTS ###
###
### Doesn't work for multi-locations (non-spatial model ok) ==> FIX THIS! (not high priority)
###
if(FALSE){
  plot.sim.all.mc(sim,"cumIncidence")
  
  # Summary across all MC iterations:
  plot.summarize.sim(sim,varname = "cumIncidence", typ="s")
  plot.summarize.sim(sim,varname = "incidence", typ="s")
  plot.summarize.sim(sim,varname = "incidence",timebucket = bucket.size, typ="s")
  
  # Summary of the simulation, agregating data in time buckets:
  sim.summary = as.data.frame(summarize.sim(sim,"cumIncidence",timebucket = NULL))
  sim.summary = as.data.frame(summarize.sim(sim,"incidence",timebucket =bucket.size))
  
  plot.Reff(sim,timebucket=7)
}


######################################################################

t1<-as.numeric(Sys.time())
msg <- paste0("\n--- SIMULATION FINISHED IN ", round( (t1-t0)/60,1), " MINUTES")
print(msg)
message(msg)
save.image(paste0("simul.RData"))
