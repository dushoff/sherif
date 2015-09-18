####################################################
####################################################
####    
###    WRAPS SPATIAL SHERIF MODEL FOR USE IN 'SYNLIK' PKG
###
###    Created 2015-09-14 by David Champredon
###
####################################################
####################################################

path.sherif.lib <- "../lib"
library(sherif,lib.loc = path.sherif.lib)
library(synlik)
library(snowfall)
source("loadParam_FCT.R")

seed.rdm <- 12345
set.seed(seed.rdm)

### Load initial model parameters
prm.model <- loadParam("param_init.csv")

### Simulation parameters for the model to be calibrated
prm.simul <- list(mc_iter = 10, # <- not used
				  horizon = 300,
				  timeStepTauLeap = 0.2,
				  calc_WIW_Re = 0,
				  seed = seed.rdm,
				  timeIdxGI = 50)

### Spatial parameters:
prm.spatial <- loadSpatialParam()
prm.model <- replicate.contact.rates(nLocations = prm.spatial[["nLocations"]], prm.model)

### Wrapping function for 'synlik'
###
sherif_spatial_wrap <- function(param, nsim, extraArgs, ...){
	
	# Unpack
	paramsModel <- extraArgs[["paramsModel"]]
	paramsSimul <- extraArgs[["paramsSimul"]]
	paramsSpatial <- extraArgs[["paramsSpatial"]]
	
	# NOTE: 'param' = model parameters to be calibrated
	# so must override model parameter given in 'paramsModel'
	for(nn in names(param)) paramsModel[nn] <- param[nn]
	
	# Overrides mc_iter by nsim
	paramsSimul[["mc_iter"]] <- nsim
	
	# Run the simulation with appropriate parameters
	sim <- rcpp_sherif_spatial(paramsSimul=paramsSimul,
							   paramsModel=paramsModel,
							   paramsSpatial=paramsSpatial)
	
	# Store the 'nsim' simulations into a list of data frames.
	# Structure of the list is:
	#
	# - location_0 ---- ts      (dataframe of time series)
	#               |__ GIbck   (dataframe of generation interval)
	#
	# - location_1 ---- ts
	#               |__ GIbck
	#
	# etc.
	
	# List of spatial locations:
	locations <- list()
	
	for(loc in 1:paramsSpatial[["nLocations"]]){
		
		# data frame of time series:
		ts <- list()
		# Generation intervals
		GIbck <- list()
		
		for(i in 1:nsim){
			df <- data.frame(time= sim[[loc]]$time[[i]],
							 cumInc= sim[[loc]]$cumIncidence[[i]],
							 deathsCum= sim[[loc]]$deathsCum[[i]])
			ts[[i]] <- df
			GIbck[[i]] <- sim[[loc]]$GIbck_gi[[i]]
		}
		locations[[loc]] <- list(ts=ts, GIbck=GIbck)
		names(locations)[loc] <- paste0("location_",loc-1)
	}
	return(locations)
}


### Statistics on the model outputs
### to inform the synthetic likelihood
###
sherif_spatial_stats <- function(x,extraArgs,...){
	
	nLocations <- length(x)
	
	for(loc in 1:nLocations){
		
		ts <- x[[loc]][["ts"]]
		GIbck <- x[[loc]][["GIbck"]]
		nsim = length(ts)
		
		# All stats
		s1 <- vector(length = nsim)
		s2 <- vector(length = nsim)
		s3 <- vector(length = nsim)
		s4 <- vector(length = nsim)
		
		for(i in 1:nsim){
			# Stat #1
			s1[i] <- mean(ts[[i]]$cumInc)	
			
			# Stat #2
			s2[i] <- sd(ts[[i]]$cumInc)
			
			# Stat #3
			s3[i] <- mean(ts[[i]]$deathsCum)
			
			# Stat #4
			s4[i] <- mean(GIbck[[i]])
		}
		
		if(loc==1) M <- cbind(s1,s2,s3,s4)
		if(loc>1) M <- cbind(M,cbind(s1,s2,s3,s4))
		
		nc <- ncol(M)
		colnames(M)[(nc-3):nc] <- paste0("loc",loc,"_",c("s1","s2","s3","s4"))
	}
	return(M)
}

### Create the synlik object
###
sherif_spatial_sl <- synlik(simulator = sherif_spatial_wrap,
							summaries = sherif_spatial_stats,
							param=c(delta = 0.61),
							extraArgs = list(
								"paramsSimul" = prm.simul, 
								"paramsModel" = prm.model,
								"paramsSpatial" = prm.spatial))

### Tests for debuging...
###
if(T){
	# A simple simulation
	zz <- simulate(sherif_spatial_sl, nsim=2)
	print(str(zz))
	ts <- zz[[1]][["ts"]]
	plot(x=ts[[1]]$time, y=ts[[1]]$cumInc, typ="l", lwd=6,
		 main="test")
	lines(x=ts[[2]]$time, y=ts[[2]]$cumInc, lwd=6)
	
	# Summary stats checks:
	ss <- sherif_spatial_stats(zz)
	print(ss)
}


stop() ##### <<<========

### Generate target data
###

delta.true <- 0.5
meanDur_latent.true <- 6

true.data <- sherif_spatial_wrap(param = c(delta = delta.true,
										   meanDur_latent = meanDur_latent.true), 
								 nsim = 1, 
								 extraArgs = list(
								 	"paramsSimul" = prm.simul, 
								 	"paramsModel" = prm.model,
								 	"paramsSpatial" = prm.spatial))

sherif_spatial_sl@data <-true.data


par(mfrow=c(3,1))

for(loc in 1:prm.spatial[["nLocations"]] ){
	true.data.ts <- true.data[[loc]][["ts"]]
	true.data.GIbck <- true.data[[loc]][["GIbck"]]
	
	plot(x=true.data.ts[[1]]$time, y=true.data.ts[[1]]$cumInc, 
		 typ="l", lwd=6,
		 main=paste0("Location #",loc, " ; Target Data - Cum Inc"))
	plot(x=true.data.ts[[1]]$time, y=true.data.ts[[1]]$deathsCum, 
		 typ="l", lwd=6,
		 main="Target Data - Cum Deaths")
	plot(density(true.data.GIbck[[1]]), 
		 typ="l", lwd=6,
		 main="Target Data - GI")
	
}


### Make these data the target one
sherif_spatial_sl@extraArgs$obsData <- sherif_spatial_sl@data

# Check normality of stats
if(T) checkNorm(sherif_spatial_sl, nsim=20)


# Range of param values explored:
# (latin hypercube style)
npts.lhs <- 20

beta_IS.rng = sample( seq(0.01, 0.7, length.out = npts.lhs) )
delta.rng = sample( seq(0.1, 0.8, length.out = npts.lhs) )
meanDur_latent.rng = sample( seq(1, 15, length.out = npts.lhs) )


t0 <- Sys.time()

sfInit(parallel = TRUE, cpu = 4)
sfLibrary(sherif,lib.loc = path.sherif.lib)
sfLibrary(synlik)

sherif_snowWrap <- function(i){
	###
	### WRAP FOR USE WITH 'SNOWFALL' PACKAGE
	###
	res <- slik(object = sherif_sl, 
				param = c(beta_IS = beta_IS.rng[i],
						  delta = delta.rng[i],
						  meanDur_latent = meanDur_latent.rng[i]), 
				multicore = F,  nsim=10)
	return(res)
}

idx.apply <- 1:npts.lhs

### Parallel execution:
sfExportAll()
system.time(xx <- sfSapply(idx.apply, sherif_snowWrap))
sfStop()




###   PLOTS   ###


df <- data.frame(beta_IS = beta_IS.rng,
				 delta = delta.rng,
				 meanDur_latent = meanDur_latent.rng,
				 synlik = xx)

par(mfrow=c(2,3))

plot(x=df$beta_IS,y=df$delta,pch=15,cex=4/log(-df$synlik))
text(x=df$beta_IS,y=df$delta,labels = round(df$synlik,0),pos=4,cex=0.7,col="grey")
abline(v=beta_IS.true,lty=2,col="red")
abline(h=delta.true,lty=2,col="red")

plot(x=df$beta_IS,y=df$meanDur_latent,pch=15,cex=4/log(-df$synlik))
text(x=df$beta_IS,y=df$meanDur_latent,labels = round(df$synlik,0),pos=4,cex=0.7,col="grey")
abline(v=beta_IS.true,lty=2,col="red")
abline(h=meanDur_latent.true,lty=2,col="red")

plot(x=df$beta_IS,y=df$synlik,pch=16)
plot(x=df$delta,y=df$synlik,pch=16)
plot(x=df$meanDur_latent,y=df$synlik,pch=16)

df[which.max(df$synlik),]

print(Sys.time()-t0)