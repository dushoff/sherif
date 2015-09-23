####################################################
####################################################
####    
###    Approximate Bayesian Computation Calibration
###
###    Created 2015-09-22 by David Champredon
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
seed <- args[2]
if(is.null(seed)) seed <- 1234
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

### Retrieve true values of parameters that generated the data
true.prm <- get.param(paramNames = names(prm2fit.val),
                      paramsModel = prm.model.true,
                      tag = "_vec")
write.csv(x = true.prm, file= paste0("true_prm.csv"),row.names = F)


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

ssfct <- select.spatial.stats(args[1])

mss <- list()
for(i in 1:npts.lhs){
  tmp <- ssfct(x = simlist[[i]],
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
prm2fit.val$ABCdist <- dist.data


### Retrieve best and top trial parameters (=shortest distances)
tmp <- rank(dist.data)
idx.best <- (tmp<=12)
idx.top <- (tmp==1)

param.best <- prm2fit.val[idx.best,]
param.best
param.top <- prm2fit.val[idx.top,]
param.top


###   save best param
filename <- paste0("best_ABC_",seed,".csv")
write.csv(x = param.best, file = filename,row.names = F)



### Retrieve true values of parameters that generated the data
true.prm <- get.param(paramNames = names(prm2fit.val),
                      paramsModel = prm.model.true,
                      tag = "_vec")


plot.pair.ABC <- function(var1,var2, true.prm, 
                          prm2fit.val,param.best, param.top){
  
  t1 <- true.prm[var1][[1]]
  t2 <- true.prm[var2][[1]]
  
  p1 <- prm2fit.val[,var1]
  p2 <- prm2fit.val[,var2]

  best1 <- param.best[,var1]
  best2 <- param.best[,var2]
  top1 <- param.top[,var1]
  top2 <- param.top[,var2]
  
  dcol <- prm2fit.val$ABCdist/max(prm2fit.val$ABCdist)
  
  plot(x=p1, y=p2, xlab=var1, ylab=var2, pch=16, cex=1, col=rgb(0,0,0,1-dcol))
  abline(v=t1,lty=2,lwd=4,col="blue")
  abline(h=t2,lty=2,lwd=4,col="blue")
  points(x=best1, y=best2, pch=1, cex=1.7, col="orange",lwd=6)
  points(x=top1, y=top2, pch=5, cex=3, lwd=4, col="red")
 
}


plot.ABC.result <- function(true.prm, 
                            prm2fit.val,
                            param.best, 
                            param.top){
  nam <- names(prm2fit.val)
  nn <- length(nam)
  q <- floor(sqrt((nn-1)*(nn-2)/2))
  par(mfrow=c(q,q))
  for(i in 1:nn)
    for(j in (i+1):nn){
      if( (!grepl("ABCdist",nam[i])) & (!grepl("ABCdist",nam[j])) ){
        plot.pair.ABC(var1=nam[i], var2=nam[j],
                      true.prm,
                      prm2fit.val,
                      param.best, 
                      param.top)
      }
    }
}

pdf("calib_result_ABC.pdf",width = 15,height = 10)
plot.ABC.result(true.prm, 
                prm2fit.val,
                param.best, 
                param.top)
dev.off()

######################################################################

t1<-as.numeric(Sys.time())
msg <- paste0("--- CALIBRATION FINISHED IN ", round( (t1-t0)/60,1), " MINUTES")
print(msg)
message(msg)
save.image(paste0("result-calib-ABC-",seed,".RData"))