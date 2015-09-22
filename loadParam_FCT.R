####
####   LOAD PARAMETERS FROM FILE INTO ENVIRONMENT
####
####   Created 2015-08-11 by David Champredon
####


loadParam <- function(filename){
	
	# unpack all model parameters from file
	# into R environment:
	prm <- read.csv(filename,header = F)
	names(prm)<-c("name","value")
	prm$name <- as.character(prm$name)
	prm.model <- list()
	for(i in 1:nrow(prm)){
		assign(as.character(prm$name[i]),prm$value[i])
		prm.model[[prm$name[i]]] <- get(prm$name[i])
	}	
	return(prm.model)
}


loadSpatialParam <- function(){
	
	### READ ALL SPATIAL PARAMETERS 
	###
	
	res <- list()
	
	# number of spatial locations:
	res[["nLocations"]] <- as.numeric(read.csv("nLocations.csv",header=FALSE))
	
	# Population size in each locations:
	res[["popLocations"]] <- as.numeric(read.csv("popLocations.csv",header=FALSE)[,1])
	
	# Distance between locations:
	M <- as.matrix(read.csv("distLocations.csv",header=FALSE))
	res[["distLocations"]] <- as.numeric(t(M))
	
	# Initial number of individuals in I, Iw and Sw:
	A <- read.csv("initLocations.csv",header=TRUE)
	res[["init_I1"]] <- A[,1]
	res[["init_Iw1"]] <- A[,2]
	res[["init_Sw1"]] <- A[,3]
	

	
	# Integrity checks
	stopifnot(res[["nLocations"]]==length(res[["popLocations"]]))
	stopifnot(res[["nLocations"]]^2==length(res[["distLocations"]]))
	stopifnot(res[["nLocations"]]==length(res[["init_I1"]]))
	stopifnot(res[["nLocations"]]==length(res[["init_Iw1"]]))
	stopifnot(res[["nLocations"]]==length(res[["init_Sw1"]]))
	
	return(res)
}


replicate.contact.rates <- function(nLocations, prm.model){
	
	### Transform single value contact rates into 
	### vector of replicated value for each location
	### (mostly used for tests)
	
	prm.model[["beta_IS"]] <- rep(prm.model[["beta_IS"]],nLocations )
	prm.model[["beta_FS"]] <- rep(prm.model[["beta_FS"]],nLocations )
	prm.model[["beta_IwS"]] <- rep(prm.model[["beta_IwS"]],nLocations )
	prm.model[["beta_ISw"]] <- rep(prm.model[["beta_ISw"]],nLocations )
	prm.model[["beta_FSw"]] <- rep(prm.model[["beta_FSw"]],nLocations )
	prm.model[["beta_IwSw"]] <- rep(prm.model[["beta_IwSw"]],nLocations )
	prm.model[["beta_HSw"]] <- rep(prm.model[["beta_HSw"]],nLocations )
	
	return(prm.model)
}


loadParamMigration <- function(filename){
	# Migration parameters:
	mig <- read.csv(filename,header=FALSE)
	return(list(migrationParams= mig[,2]))
}
