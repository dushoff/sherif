####################################################
####################################################
####    
###    RUNS SHERIF MODEL WITHIN R (not command line c++)
###
###    Created 2015-08-19 by David Champredon
###
####################################################
####################################################

library(plyr)
library(reshape2)
library(snowfall)

### R package for SHERIF is installed locally:
###
path.sherif.lib <- "./lib"

library(sherif, lib.loc = path.sherif.lib)

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

