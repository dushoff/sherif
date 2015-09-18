####################################################
####################################################
####    
###    WRAPS SHERIF MODEL FOR USE IN 'SYNLIK' PKG
###
###    Created 2015-08-11 by David Champredon
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


### Wrapping function for 'synlik'
###
sherif_wrap <- function(param, nsim, extraArgs, ...){
	
	# Unpack
	paramsModel <- extraArgs[["paramsModel"]]
	paramsSimul <- extraArgs[["paramsSimul"]]
	
	# NOTE: 'param' = model parameters to be calibrated
	# so must override model parameter given in 'paramsModel'
	for(nn in names(param)) paramsModel[nn] <- param[nn]
	
	# Overrides mc_iter by nsim
	paramsSimul[["mc_iter"]] <- nsim
	
	# Run the simulation with appropriate parameters
	sim <- rcpp_sherif(paramsSimul=paramsSimul,
					   paramsModel=paramsModel)
	
	# Store the 'nsim' simulations into a list of data frames:
	
	# data frame of time series:
	ts <- list()
	# Generation intervals
	GIbck <- list()
	
	for(i in 1:nsim){
		df <- data.frame(time= sim$time[[i]],
						 cumInc= sim$cumIncidence[[i]],
						 deathsCum= sim$deathsCum[[i]])
		ts[[i]] <- df
		GIbck[[i]] <- sim$GIbck_gi[[i]]
	}
	return(list(ts=ts, GIbck=GIbck))
}

### Statistics on the model outputs
### to inform the synthetic likelihood
###
sherif_stats <- function(x,extraArgs,...){
	
	ts <- x[["ts"]]
	GIbck <- x[["GIbck"]]
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
	return(cbind(s1,s2,s3,s4))
}

### Create the synlik object
###
sherif_sl <- synlik(simulator = sherif_wrap,
					summaries = sherif_stats,
					param=c(beta_IS = 0.33),
					extraArgs = list(
						"paramsSimul" = prm.simul, 
						"paramsModel" = prm.model))

if(FALSE){
	zz <- simulate(sherif_sl, nsim=1)
	ts <- zz[["ts"]]
	plot(x=ts[[1]]$time, y=ts[[1]]$cumInc, typ="l", lwd=6,
		 main="test")
}

### Generate target data
###

beta_IS.true <- 0.20
delta.true <- 0.5
meanDur_latent.true <- 6

true.data <- sherif_wrap(param = c(beta_IS = beta_IS.true,
								   delta = delta.true,
								   meanDur_latent = meanDur_latent.true), 
						 nsim = 1, 
						 extraArgs = list(
						 	"paramsSimul" = prm.simul, 
						 	"paramsModel" = prm.model))

sherif_sl@data <-true.data

true.data.ts <- true.data[["ts"]]
true.data.GIbck <- true.data[["GIbck"]]

###   MCMC inference
###
if(FALSE){
	mcmc_sherif_sl <- smcmc(sherif_sl, 
							initPar = c(0.2, 0.5, 6),
							niter = 20, 
							burn = 7,
							priorFun = function(input, ...) sum(input), 
							propCov = diag(c(0.1, 0.1, 0.1))^2, 
							nsim = 10)
	
	plot(mcmc_sherif_sl@chains[,3])
}



par(mfrow=c(3,1))

plot(x=true.data.ts[[1]]$time, y=true.data.ts[[1]]$cumInc, 
	 typ="l", lwd=6,
	 main="Target Data - Cum Inc")
plot(x=true.data.ts[[1]]$time, y=true.data.ts[[1]]$deathsCum, 
	 typ="l", lwd=6,
	 main="Target Data - Cum Deaths")
plot(density(true.data.GIbck[[1]]), 
	 typ="l", lwd=6,
	 main="Target Data - GI")

### Make these data the target one
sherif_sl@extraArgs$obsData <- sherif_sl@data

# Check normality of stats
if(F) checkNorm(sherif_sl, nsim=20)


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