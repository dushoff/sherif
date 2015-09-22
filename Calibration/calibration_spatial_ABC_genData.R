####################################################
####################################################
####    
###    WRAPS SHERIF MODEL FOR USE IN 'SYNLIK' PKG
###
###    Created 2015-08-11 by David Champredon
###
####################################################
####################################################

t0 <- as.numeric(Sys.time())

### Command line arguments
### args[1]: summary statistics function to be used (integer=1,2,3,4)
### args[2]: seed for random number generator
###

# args <- commandArgs(trailingOnly = TRUE)

### R package for SHERIF is installed locally:
path.sherif.lib <- "../lib"

library(sherif,lib.loc = path.sherif.lib)
library(snowfall)
# detect the number of CPU cores
library(parallel) 
n.cpu <- detectCores(all.tests = TRUE, logical = FALSE)   
source("loadParam_FCT.R")
source("wrap_FCT.R")
source("synlik_stats.R")
source("inference_FCT.R")

### Load initial model parameters
prm.model <- loadParam("param_init.csv")

### Simulation parameters for the model to be calibrated:
prm.simul <- loadParam("param_simul.csv")
set.seed(prm.simul[["seed"]])
horizon <- prm.simul[["horizon"]]

### Spatial parameters:
prm.spatial <- loadSpatialParam()
prm.model <- replicate.contact.rates(nLocations = prm.spatial[["nLocations"]], prm.model)
prm.model <- c(prm.model,loadParamMigration("gravity_cst.csv"))

### Generate target data
### (parameter values in 'param_true4fit.csv')
###
seed <- 1234
file.data <- paste0("generated-target-data.pdf")
if(!is.na(seed)) file.data <- paste0("generated-target-data_",seed,".pdf")
pdf(file.data, width = 12, height=10)
source("generate-target-data-spatial_ABC.R")
dev.off()



###########################################
###
###    CALIBRATION
###
###########################################

### Parameters for the ABC trials:
lhsp <- read.csv("ABCParam.csv")
npts.lhs <- lhsp$value[lhsp$var=="nlhs"]
nmc <- lhsp$value[lhsp$var=="nmc"]
message(paste0("n.lhs=",npts.lhs," ; nmc=",nmc,
               " ==> Total of ", nmc*npts.lhs, " iterations"))

### Latin hypercube sampling of desired parameters.
### (just add the name of the variable in the data frame
###  to include it to the list of parameters to fit)
###
prm2fit.val <- data.frame(beta_IS_vec1 = sample( seq(0.01, 0.7,length.out = npts.lhs)),
                          beta_IS_vec2 = sample( seq(0.01, 0.7,length.out = npts.lhs)),
                          beta_IS_vec3 = sample( seq(0.01, 0.7,length.out = npts.lhs)),
                          #migrationParams_vec2 = sample( seq(2, 7,length.out = npts.lhs)),
                          meanDur_latent = sample( seq(1, 15,length.out = npts.lhs))
)
prior.range <- rbind(apply(prm2fit.val,MARGIN = 2, FUN = min),
                     apply(prm2fit.val,MARGIN = 2, FUN = max))

### Retrieve true values of parameters that generated the data
true.prm <- get.param(paramNames = names(prm2fit.val),
                      paramsModel = prm.model.true,
                      tag = "_vec")

# Summary statistics of the (generated) data
SS.data <- sherif_spatial_stats_1(x = true.data2,
                                  extraArgs = list(horizon=horizon))

all.baseline.prm <-  list(paramsSimul = prm.simul,
                          paramsModel = prm.model,
                          paramsSpatial = prm.spatial)

### Calculate summary stats for every values
### of the LHS sampling
###
ABC_snowWrap <- function(i){
  ### WRAP FOR USE WITH 'SNOWFALL' PACKAGE
  ###
  newprm <- as.numeric(prm2fit.val[i,])
  names(newprm) <- colnames(prm2fit.val)
  res <- sherif_spatial_wrap(param = newprm,
                                      extraArgs = all.baseline.prm,
                                      nsim = nmc)
  return(res)
}

idx.apply <- 1:npts.lhs

### Parallel execution:
# Snowfall setup for parallel computing:
sfInit(parallel = TRUE, cpu = n.cpu)
sfLibrary(sherif, lib.loc = path.sherif.lib)
sfExportAll()
system.time(simlist <- sfSapply(idx.apply, ABC_snowWrap,simplify = F))
sfStop()

### Mean of the summary stats across all Monte Carlo samples
mss <- list()
for(i in 1:npts.lhs){
  tmp <- sherif_spatial_stats_1(x = simlist[[i]],
                                extraArgs = list(horizon=horizon))
  mss[[i]] <- apply(tmp,MARGIN = 2,FUN = mean)
}

### Distance function
dist.relative <- function(x,y,weights){
  # Relative distance element-wise, with weights
  stopifnot(length(x)==length(y))
  stopifnot(length(x)==length(weights))
  res <- NA
  if(!(all(y)==0)) res <- sqrt(sum(weights*(x/y-1)^2))
  return(res)
}

### Distance of the vectors
### of summary stats mean calculated with trial
### parameters to the vector of summary stats  
### from the data to fit
dist.data <- numeric(npts.lhs)
weigths <- rep(1.0,length(SS.data))
for(i in 1:npts.lhs){
  dist.data[i] <-dist.relative(mss[[i]],as.numeric(SS.data),weigths) 
}


### Retrieve best and top trial parameters (=shortest distances)
tmp <- order(dist.data)
idx.best <- (tmp<=3)
idx.top <- (tmp==1)

param.best <- prm2fit.val[idx.best,]
param.best$ABCdist <- dist.data[idx.best]
param.best

param.top <- prm2fit.val[idx.top,]
param.top



### Retrieve true values of parameters that generated the data
true.prm <- get.param(paramNames = names(prm2fit.val),
                      paramsModel = prm.model.true,
                      tag = "_vec")


plot.pair.ABC <- function(var1,var2, true.prm){
  var1="beta_IS_vec1"
  var2 = "meanDur_latent"
  t1 <- true.prm[var1]
  t1 <- true.prm[var1]
}

plot.ABC.result <- function(param.best, param.top,true.prm){
  
}

######################################################################

t1<-as.numeric(Sys.time())
msg <- paste0("--- CALIBRATION FINISHED IN ", round( (t1-t0)/60,1), " MINUTES")
print(msg)
message(msg)
# save.image(paste0("result-calib-",synlik_stats_type,".RData"))