####################################################
####################################################
####    
###    WRAPS SHERIF MODEL FOR USE IN 'SYNLIK' PKG
###
###    Created 2015-08-11 by David Champredon
###
####################################################
####################################################


### Command line arguments
### args[1]: summary statistics function to be used (integer=1,2,3,4)
### args[2]: seed for random number generator
###
args <- commandArgs(trailingOnly = TRUE)


t0 <- as.numeric(Sys.time())

### R package for SHERIF is installed locally:
path.sherif.lib <- "../lib"

library(sherif,lib.loc = path.sherif.lib)
library(synlik)
library(snowfall)
library(methods)
library(parallel) # <-- just to detect the number of CPU cores

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


### Select the type of summary statistics for synlik here:
### (functions defined in 'synlik_stats.R')
###
synlik_stats_type <- args[1]
sherif_stats <- select.synlik.stats(synlik_stats_type)


### Create the synlik object
###
sherif_sl <- synlik(simulator = sherif_wrap,
                    summaries = sherif_stats,
                    param=c(beta_IS = 0.33),
                    extraArgs = list(
                      "paramsSimul" = prm.simul, 
                      "paramsModel" = prm.model))


### Generate target data
### (parameter values in 'param_true4fit.csv')
###
file.data <- paste0("generated-target-data.pdf")
if(!is.na(seed)) file.data <- paste0("generated-target-data_",seed,".pdf")

pdf(file.data, width = 12, height=10)
source("generate-target-data.R")

# Update the synlik object with
# the new 'horizon' (latest date of observed data)
# which is needed for summary stats
sherif_sl <- synlik(simulator = sherif_wrap,
                    summaries = sherif_stats,
                    param=c(beta_IS = 0.33),
                    extraArgs = list(
                      "paramsSimul" = prm.simul, 
                      "paramsModel" = prm.model,
                      "horizon" = horizon))

sherif_sl@data <- true.data2
sherif_sl@extraArgs$obsData <- sherif_sl@data


###########################################
###
###    CALIBRATION
###
###########################################



n.cpu <- detectCores(all.tests = TRUE, logical = FALSE)   # as.integer(read.table("ncpu.csv"))
lhsp <- read.csv("lhsParam.csv")
npts.lhs <- lhsp$value[lhsp$var=="nlhs"]
nsim <- lhsp$value[lhsp$var=="nsim"]

# Check normality of stats
if(TRUE) {
  pdf(paste0("check_norm_",synlik_stats_type,".pdf"))
  checkNorm(sherif_sl, nsim=nsim)
  dev.off()
}

message(paste0("n.lhs=",npts.lhs," ; nsim=",nsim,
               " ==> total of ", nsim*npts.lhs, " iterations"))

### Latin hypercube sampling of desired parameters.
### (just add the name of the variable in the data frame
###  to include it to the list of parameters to fit)
###
prm2fit.val <- data.frame(beta_IS = sample( seq(0.01, 0.7,length.out = npts.lhs)),
                          beta_FS = sample( seq(0.01, 0.7,length.out = npts.lhs)),
                          meanDur_funeral =  sample( seq(0.5,4,length.out = npts.lhs)),
                          #delta = sample( seq(0.1, 0.9,length.out = npts.lhs)),
                          meanDur_latent = sample( seq(1, 15,length.out = npts.lhs))
)
prior.range <- rbind(apply(prm2fit.val,MARGIN = 2, FUN = min),
                     apply(prm2fit.val,MARGIN = 2, FUN = max))

### Retrieve true values of parameters that generated the data
true.prm <- prm.model.true[names(prm2fit.val)]


### Calibration result with Latin Hypercube Sampling:
###
file.plot <- paste0("calib_results_",synlik_stats_type,".pdf")
if(!is.na(seed)) file.plot <- paste0("calib_results_",synlik_stats_type,"_",seed,".pdf")

pdf(file.plot, width=15, height=12)

calib.lhs <- inference_LHS(prm2fit.val,
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
# write.table(post, sep = ",",
#             file=file.post.table,
#             quote = FALSE, col.names = TRUE)
write.csv(post, file=file.post.table)


### Save in file max synlik 
### (for further comparison when several summary stats are run)
###
file.maxlik <- paste0("calib_results_maxLik_",synlik_stats_type,".csv")
if(!is.na(seed)) file.maxlik <- paste0("calib_results_maxLik_",synlik_stats_type,"_",seed,".csv")

write.table(diff.from.true(true.prm,max.synlik), sep = ",",
            file=file.maxlik,
            quote = FALSE, col.names = FALSE)




#### WARNING : NO PARALLEL COMPUTATION IMPLEMENTED FOR MCMC (FIX THIS)
#### ==> CAN BE LONG!
if(FALSE){
  res.mcmc <- inference_mcmc(synlik.object = sherif_sl,
                             mcmc.iter = 9, 
                             prm.init = c(0.4, 0.3, 9),
                             mcmc.warmup = 2,
                             nsim = 3,
                             do.plot = TRUE)
}


######################################################################

t1<-as.numeric(Sys.time())
msg <- paste0("--- CALIBRATION FINISHED IN ", round( (t1-t0)/60,1), " MINUTES")
print(msg)
message(msg)
dev.off()
save.image(paste0("result-calib-",synlik_stats_type,".RData"))