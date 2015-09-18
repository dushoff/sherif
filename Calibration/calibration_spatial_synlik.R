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

args <- commandArgs(trailingOnly = TRUE)

### R package for SHERIF is installed locally:
path.sherif.lib <- "../Rcpp-sherif/library"

library(sherif,lib.loc = path.sherif.lib)
library(synlik)
library(snowfall)
library(methods)
# detect the number of CPU cores
library(parallel) 
n.cpu <- detectCores(all.tests = TRUE, logical = FALSE)   
source("loadParam_FCT.R")
source("wrap_FCT.R")
source("synlik_stats.R")
source("inference_FCT.R")

seed <- args[2]
tmp <- ifelse(is.na(seed),1234,seed) 
print(paste("SEED=",tmp))
set.seed(tmp)

### Load initial model parameters
prm.model <- loadParam("param_init.csv")

### Simulation parameters for the model to be calibrated:
prm.simul <- loadParam("param_simul.csv")
set.seed(prm.simul[["seed"]])

### Spatial parameters:
prm.spatial <- loadSpatialParam()
prm.model <- replicate.contact.rates(nLocations = prm.spatial[["nLocations"]], prm.model)
prm.model <- c(prm.model,loadParamMigration("gravity_cst.csv"))

### Select the type of summary statistics for synlik here:
### (functions defined in 'synlik_stats.R')
###
synlik_stats_type <- args[1]
sherif_spatial_stats <- select.synlik.spatial.stats(synlik_stats_type)

### Create the synlik object
###
sherif_spatial_sl <- synlik(simulator = sherif_spatial_wrap,
                            summaries = sherif_spatial_stats,
                            param=c(delta = 0.33),
                            extraArgs = list(
                              "paramsSimul" = prm.simul, 
                              "paramsModel" = prm.model,
                              "paramsSpatial" = prm.spatial))

### Generate target data
### (parameter values in 'param_true4fit.csv')
###
file.data <- paste0("generated-target-data.pdf")
if(!is.na(seed)) file.data <- paste0("generated-target-data_",seed,".pdf")
pdf(file.data, width = 12, height=10)
source("generate-target-data-spatial.R")

# Update the synlik object with
# the new 'horizon' (latest date of observed data)
# which is needed for summary stats
sherif_spatial_sl <- synlik(simulator = sherif_spatial_wrap,
                            summaries = sherif_spatial_stats,
                            param=c(delta = 0.33),
                            extraArgs = list(
                              "paramsSimul" = prm.simul, 
                              "paramsModel" = prm.model,
                              "paramsSpatial" = prm.spatial,
                              "horizon" = horizon)) # <-- 'horizon' added
sherif_spatial_sl@data <- true.data2
sherif_spatial_sl@extraArgs$obsData <- sherif_spatial_sl@data

### Check if summary stats are OK:
if(TRUE){
  xx <- sherif_spatial_wrap(param = c(delta = prm.model.true[["delta"]]), # <-- does not matter which param chosen here. Just to comply with function required signature
                            nsim = 5, 
                            extraArgs = list(
                              "paramsSimul" = prm.simul.true, 
                              "paramsModel" = prm.model.true,
                              "paramsSpatial" = prm.spatial.true,
                              "horizon" = horizon),
                            seed = seed)
  ss <- sherif_spatial_stats_1(x=xx,extraArgs = list(
    "paramsSimul" = prm.simul.true, 
    "paramsModel" = prm.model.true,
    "paramsSpatial" = prm.spatial.true,
    "horizon" = horizon))
  print(ss)
  # stop if there are NAs --> must fix the summary stats:
  stopifnot(sum(is.na(ss))==0)
}

###########################################
###
###    CALIBRATION
###
###########################################


lhsp <- read.csv("lhsParam.csv")
npts.lhs <- lhsp$value[lhsp$var=="nlhs"]
nsim <- lhsp$value[lhsp$var=="nsim"]

# Check normality of stats
if(TRUE) {
  pdf(paste0("check_norm_",synlik_stats_type,".pdf"))
  checkNorm(sherif_spatial_sl, nsim=nsim)
  dev.off()
}

message(paste0("n.lhs=",npts.lhs," ; nsim=",nsim,
               " ==> total of ", nsim*npts.lhs, " iterations"))

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


### Calibration result with Latin Hypercube Sampling:
###
file.plot <- paste0("calib_results_",synlik_stats_type,".pdf")
if(!is.na(seed)) file.plot <- paste0("calib_results_",synlik_stats_type,"_",seed,".pdf")

pdf(file.plot, width=15, height=12)

calib.lhs <- inference_spatial_LHS(prm2fit.val,
                                   n.cpu = n.cpu, nsim = nsim,
                                   path.sherif.lib = path.sherif.lib,
                                   do.plot = TRUE, 
                                   true.prm = true.prm)

max.synlik = calib.lhs[["bestfit"]]
df = calib.lhs[["alltrials"]]

### Posteriors
###
file.post <- paste0("calib_results_post_",synlik_stats_type,".pdf")
if(!is.na(seed)) file.post <- paste0("calib_results_post_",synlik_stats_type,"_",seed,".pdf")

pdf(file.post, width=12,height=12)
post <- posteriors(df,true.prm, do.plot=TRUE,prior.range = prior.range)
dev.off()

file.post.table <- paste0("calib_results_post_",synlik_stats_type,".csv")
if(!is.na(seed)) file.post.table <- paste0("calib_results_post_",synlik_stats_type,"_",seed,".csv")
write.csv(post,file=file.post.table)


### Save in file max synlik 
### (for further comparison when several summary stats are run)
###
file.maxlik <- paste0("calib_results_maxLik_",synlik_stats_type,".csv")
if(!is.na(seed)) file.maxlik <- paste0("calib_results_maxLik_",synlik_stats_type,"_",seed,".csv")

write.table(diff.from.true(true.prm,max.synlik), sep = ",",
            file=file.maxlik,
            quote = FALSE, col.names = FALSE)




######################################################################

t1<-as.numeric(Sys.time())
msg <- paste0("--- CALIBRATION FINISHED IN ", round( (t1-t0)/60,1), " MINUTES")
print(msg)
message(msg)
dev.off()
save.image(paste0("result-calib-",synlik_stats_type,".RData"))