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
prm.model.true <- replicate.contact.rates(nLocations = prm.spatial.true[["nLocations"]], prm.model.true)



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
nLocations = prm.spatial.true[["nLocations"]]
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
	
	true.data2[[loc]] <- list(ts = list(true.data[[loc]][["ts"]][[1]][t.idx,]), # <-- the "list" is need to keep the original object structure
							  GIbck = true.data[[loc]]$GIbck)
}


### Make these data the target one
###
sherif_spatial_sl@data <- true.data2 # true.data2
sherif_spatial_sl@extraArgs$obsData <- sherif_spatial_sl@data


### Plot the generated target data:
###
ts.plot <- function(ts,var, title){
	head(ts)
	typ <- "s"
	if(length(ts$time)>50) typ <- "l"
	plot(x=ts$time, y=ts[,var], typ=typ, lwd=6, 
		 main=title, xlab="time", ylab=var)
	grid()
}

if(plot.true.data)
{
	nLocations <- length(true.data2)
	
	for(loc in 1:nLocations){
		ts.true <- true.data2[[loc]][["ts"]][[1]]
		gi.true <- true.data2[[loc]][["GIbck"]][[1]]
		
		par(mfrow=c(3,2))
		title <- paste0("Location #",loc," Target Data - ")
		ts.plot(ts.true,"cumInc", paste(title,"Cum Inc"))
		ts.plot(ts.true,"cumInc_hcw", paste(title,"Cum Inc HCW"))
		ts.plot(ts.true,"deathsCum", paste(title,"Cum Deaths"))
		ts.plot(ts.true,"deathsCum_hcw", paste(title,"Cum Deaths HCW"))
		hist(gi.true, 
			 breaks = seq(0,max(ceiling(gi.true))+4,by=1),
			 col = "lightgrey",
			 main = paste(title,"Histogram bckwd GI"))
		
	}
	dev.off()
}
