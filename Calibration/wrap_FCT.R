
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
	
	# Overrides the seed for random number generator
	if(!is.null(seed)) paramsSimul[["seed"]] <- as.integer(seed)
	
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
						 cumInc_hcw= sim$cumIncidence_hcw[[i]],
						 deathsCum= sim$deathsCum[[i]],
						 deathsCum_hcw= sim$deathsCum_hcw[[i]])
		ts[[i]] <- df
		GIbck[[i]] <- sim$GIbck_gi[[i]]
	}
	return(list(ts=ts, GIbck=GIbck))
}


update.param <- function(param, paramsModel, tag){
	
	### 'param' = model parameters to be calibrated
	### so must override model parameter given in 'paramsModel'
	###
	### But need special treatment for parameters in vector format
	### (e.g. contact rates for each location)
	
	# identify parameters that are element of a vector
	# (for now, only spatial related, so look for "loc" string)
	tag <- "_loc"
	idxvec <- grepl(pattern = tag, x = names(param))
	
	# overides scalar (non-vector) parameter values 
	scalarname <- names(param)[!idxvec]
	for(nn in scalarname) paramsModel[nn] <- param[nn]
	
	if(sum(idxvec)>0){
		# deal with vector elements (a bit more fidly)
		vecname <- 	names(param)[idxvec]
		pos <- as.integer(gregexpr(pattern = "_loc",vecname))
		vecname2 <- substr(vecname,start=1,stop = pos-1)
		# element number in the vector
		elemnum <- as.integer(substr(vecname,start=pos+nchar(tag),stop = pos+nchar(tag)))
		
		# Replace the parameter value of the vector element:
		for(i in 1:length(vecname2)){
			vn <- vecname2[i]
			paramsModel[vn][[1]][elemnum[i]] <- param[vecname[i]]
		}
	}
	return(paramsModel)
}


get.param <- function(paramNames, paramsModel, tag){
	
	
	res <- vector()
	k = 1
	tag <- "_loc"
	idxvec <- grepl(pattern = tag, x = paramNames)
	
	# overides scalar (non-vector) parameter values 
	scalarname <- paramNames[!idxvec]
	for(nn in scalarname) {
		res[k] <- paramsModel[nn]
		names(res)[k] <- nn
		k = k+1
	}
	
	if(sum(idxvec)>0){
		# deal with vector elements (a bit more fidly)
		vecname <- 	paramNames[idxvec]
		pos <- as.integer(gregexpr(pattern = "_loc",vecname))
		vecname2 <- substr(vecname,start=1,stop = pos-1)
		# element number in the vector
		elemnum <- as.integer(substr(vecname,start=pos+nchar(tag),stop = pos+nchar(tag)))
		
		# Replace the parameter value of the vector element:
		for(i in 1:length(vecname2)){
			vn <- vecname2[i]
			res[k] <- paramsModel[vn][[1]][elemnum[i]] 
			names(res)[k] <- vecname[i]
			k <- k+1
		}
	}
	return(res)
}




### Wrapping function for the SPATIAL model
###
sherif_spatial_wrap <- function(param, nsim, extraArgs, ...){
	
	# Unpack
	paramsModel <- extraArgs[["paramsModel"]]
	paramsSimul <- extraArgs[["paramsSimul"]]
	paramsSpatial <- extraArgs[["paramsSpatial"]]
	
	# 	### 'param' = model parameters to be calibrated
	# 	### so must override model parameter given in 'paramsModel'
	# 	###
	# 	### But need special treatment for parameters in vector format
	# 	### (e.g. contact rates for each location)
	# 	
	# 	# identify parameters that are element of a vector
	# 	# (for now, only spatial related, so look for "loc" string)
	# 	tag <- "_loc"
	# 	idxvec <- grepl(pattern = tag, x = names(param))
	# 	
	# 	# overides scalar (non-vector) parameter values 
	# 	scalarname <- names(param)[!idxvec]
	# 	for(nn in scalarname) paramsModel[nn] <- param[nn]
	# 	
	# 	# deal with vector elements (a bit more fidly)
	# 	vecname <- 	names(param)[idxvec]
	# 	pos <- as.integer(gregexpr(pattern = "_loc",vecname))
	# 	vecname2 <- substr(vecname,start=1,stop = pos-1)
	# 	# element number in the vector
	# 	elemnum <- as.integer(substr(vecname,start=pos+nchar(tag),stop = pos+nchar(tag)))
	# 	
	# 	# Replace the parameter value of the vector element:
	# 	for(i in 1:length(vecname2)){
	# 		vn <- vecname2[i]
	# 		paramsModel[vn][[1]][elemnum[i]] <- param[vecname[i]]
	# 	}
	
	paramsModel <- update.param(param = param,
								paramsModel = paramsModel,
								tag="_loc")
	
	
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
							 deathsCum= sim[[loc]]$deathsCum[[i]],
							 cumInc_hcw= sim[[loc]]$cumIncidence_hcw[[i]],
							 deathsCum_hcw= sim[[loc]]$deathsCum_hcw[[i]])
			ts[[i]] <- df
			GIbck[[i]] <- sim[[loc]]$GIbck_gi[[i]]
		}
		locations[[loc]] <- list(ts=ts, GIbck=GIbck)
		names(locations)[loc] <- paste0("location_",loc-1)
	}
	return(locations)
}


