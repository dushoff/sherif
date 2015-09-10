##############################
###
###   Generate target data
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

### simulate one iteration to generate target data
###
true.data <- sherif_wrap(param = c(beta_IS = prm.model.true[["beta_IS"]]), # <-- does not matter which param chosen here. Just to comply with function required signature
                         nsim = 1, 
                         extraArgs = list(
                           "paramsSimul" = prm.simul.true, 
                           "paramsModel" = prm.model.true,
                           "horizon" = horizon),
                         seed = seed)


### Truncate data length (to emulate early epidemic)
### and thin the observations
###
obs.horizon <- prm.gen.data$value[prm.gen.data$name=="horizon"]
n.obs <- prm.gen.data$value[prm.gen.data$name=="n.obs"]
tmpdf <- true.data[["ts"]][[1]]
tmp = subset(tmpdf, time <= obs.horizon) 
t.idx <- sort(sample(1:length(tmp$time),size = n.obs,replace = F))
true.data2 <- list(ts = list(true.data[["ts"]][[1]][t.idx,]), # <-- the "list" is need to keep the original object structure
                   GIbck = true.data$GIbck)


### Make these data the target one
###
sherif_sl@data <- true.data2 # true.data2
sherif_sl@extraArgs$obsData <- sherif_sl@data
# sherif_sl.tmp@data <- true.data2 # true.data2
# sherif_sl.tmp@extraArgs$obsData <- sherif_sl.tmp@data


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
  ts.true <- true.data2[["ts"]][[1]]
  gi.true <- true.data2[["GIbck"]][[1]]
  
  par(mfrow=c(3,2))
  ts.plot(ts.true,"cumInc", "Target Data - Cum Inc")
  ts.plot(ts.true,"cumInc_hcw", "Target Data - Cum Inc HCW")
  ts.plot(ts.true,"deathsCum", "Target Data - Cum Deaths")
  ts.plot(ts.true,"deathsCum_hcw", "Target Data - Cum Deaths HCW")
  hist(gi.true, 
       breaks = seq(0,max(ceiling(gi.true))+4,by=1),
       col = "lightgrey",
       main = "Target Data - Histogram bckwd GI")
  dev.off()
}
