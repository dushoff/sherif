####################################################
####################################################
####    
###    RUNS SHERIF MODEL WITHIN R (not command line c++)
###
###    Created 2015-08-19 by David Champredon
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
source("run_sherif_FCT.R")
source("./Calibration/loadParam_FCT.R")

try(system("mkdir NIH_format",intern = F),silent = TRUE)

### Load simulation parameters
###
prm.simul <- loadParam("param_simul.csv")
prm.model <- loadParam("param_model.csv")


### Run the SHERIF simulations
###
ncpus <- 4
sim <- run.sherif.parallel(prm.simul, prm.model, ncpus, path.sherif.lib)


### Save output in NIH's format
proba.inc <- seq(0,0.12, by=0.001)
proba.cuminc <- seq(0,1, by=0.01)
prediction.date <- read.table("prediction_date.csv")[1,1]
bucket.size <- 7
filesuffix <- "REG_1"

NIH.format.sherif.output(sim, 
                         popsize = prm.model$popSize,
                         proba.inc,
                         proba.cuminc,
                         prediction.date,
                         bucket.size,
                         filesuffix)

### PLOTS ###
if(TRUE){
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
