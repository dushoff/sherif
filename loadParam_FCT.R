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
