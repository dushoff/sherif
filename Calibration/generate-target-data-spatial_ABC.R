##############################
###
###   Generate target data
###   for SHERIF spatial model
###
##############################

plot.true.data <- TRUE


### Load the 'true' parameters that will
### generate the target data 
###
prm.model.true <- loadParam("param_true4fit.csv")
prm.simul.true <- loadParam("param_simul_true4fit.csv")

### Load parameters to thin simulated data
prm.gen.data <- read.csv("generate-target-data-param.csv")
horizon <- prm.gen.data$value[prm.gen.data$name=="horizon"]

### Load spatial parameters:
prm.spatial.true <- loadSpatialParam()
nLocations = prm.spatial.true[["nLocations"]]
prm.model.true <- replicate.contact.rates(nLocations = nLocations, prm.model.true)
prm.model.true <- c(prm.model.true,loadParamMigration("gravity_cst.csv"))

### simulate one iteration to generate target data
###
true.data <- sherif_spatial_wrap(param = c(delta = prm.model.true[["delta"]]), # <-- does not matter which param chosen here. Just to comply with function required signature
								 nsim = 1, 
								 extraArgs = list(
								 	"paramsSimul" = prm.simul.true, 
								 	"paramsModel" = prm.model.true,
								 	"paramsSpatial" = prm.spatial.true,
								 	"horizon" = horizon),
								 seed = seed)


### Truncate data length (to emulate early epidemic)
### and thin the observations
###
obs.horizon <- prm.gen.data$value[prm.gen.data$name=="horizon"]
n.obs <- prm.gen.data$value[prm.gen.data$name=="n.obs"]

true.data2 <- list()

for(loc in 1:nLocations){
	# Retrieve the time series of generated data:
	tmpdf <- true.data[[loc]][["ts"]][[1]]
	# Truncate before observation horizon:
	tmp = subset(tmpdf, time <= obs.horizon) 
	# Draw random observation dates between now and observation horizon:
	t.idx <- sort(sample(1:length(tmp$time),size = n.obs,replace = F))
	
	true.data2[[loc]] <- list(ts = list(true.data[[loc]][["ts"]][[1]][t.idx,]), # <-- the "list" is needed to keep the original object structure
							  GIbck = true.data[[loc]]$GIbck,
							  time_firstCase = true.data[[loc]]$time_firstCase)
}



### Plot the generated target data:
###
ts.plot <- function(ts,var, title, time.firstcase=NULL){
	head(ts)
	typ <- "s"
	if(length(ts$time)>50) typ <- "l"
	plot(x=ts$time, y=ts[,var], typ=typ, lwd=6, 
		 xlim=c(0,max(ts$time)),
		 main=title, xlab="time", ylab=var)
	if(!is.null(time.firstcase)){
		abline(v=time.firstcase,lty=2)
		text(x=time.firstcase,
			 y=max(ts[,var])*0.9,
			 labels = paste("Time first case =",time.firstcase),
			 pos = 4)
	}
	grid()
}

if(plot.true.data)
{
	nLocations <- length(true.data2)
	
	for(loc in 1:nLocations){
		ts.true <- true.data2[[loc]][["ts"]][[1]]
		gi.true <- true.data2[[loc]][["GIbck"]][[1]]
		t.1stCase<- true.data2[[loc]][["time_firstCase"]][[1]]
		
		par(mfrow=c(3,2))
		title <- paste0("Location #",loc," Target Data - ")
		ts.plot(ts.true,"cumInc", paste(title,"Cum Inc"),time.firstcase=t.1stCase)
		ts.plot(ts.true,"cumInc_hcw", paste(title,"Cum Inc HCW"))
		ts.plot(ts.true,"deathsCum", paste(title,"Cum Deaths"))
		ts.plot(ts.true,"deathsCum_hcw", paste(title,"Cum Deaths HCW"))
		try(hist(gi.true, 
			 breaks = seq(0,max(ceiling(gi.true))+4,by=1),
			 col = "lightgrey",
			 main = paste(title,"Histogram bckwd GI")),silent = TRUE)
		
	}
}
