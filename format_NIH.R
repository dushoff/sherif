####################################################################
###
###    FORMATTING OF RESULTS FOR NIH'S RAPIDD EBOLA CHALLENGE
###
###    Created 2015-08-21 by David Champredon
###
####################################################################




save.NIH.format<-function(x,popsize,filename,unscale=TRUE){
  df <- x[["proba.interval"]]
  if(unscale) df$interval <- df$interval/popsize
  write.csv(df,file = filename)
}

NIH.format.sherif.output <- function(
	sim, popsize,
	proba.inc, proba.cuminc, prediction.date,
	bucket.size,filesuffix, folder
){
  par(mfrow=c(2,2))
  
  ### Rescale to indiviual units
  int.proba.inc <- round(proba.inc*popsize)
  int.proba.cuminc <- round(proba.cuminc*popsize)
  
  
  ### Incidence at a given week
  inc <- distrib.output.sherif(sim, varname="incidence",
                               interv.proba = int.proba.inc,
                               t=prediction.date, timebucket = bucket.size,
                               do.plot=TRUE)
  save.NIH.format(inc,popsize = popsize, 
                  filename = paste0(folder,"inc_proba_",filesuffix,".csv"))
  
  
  ### Final size
  last.date <- floor(max(sim$time[[1]]))
  cumInc <- distrib.output.sherif (sim, varname="cumIncidence",
                                   interv.proba = int.proba.cuminc,
                                   t=last.date,  timebucket = NULL,
                                   do.plot=TRUE)
  save.NIH.format(cumInc,popsize = popsize, 
                  filename = paste0(folder,"finalsize_proba_",filesuffix,".csv"))
  
  
  ### Timing peak incidence
  tmax <- time.peakIncidence(sim, timebucket = bucket.size)/bucket.size #  <-- 7 because weekly forcast
  horizon <- max(tmax)+1
  h <- hist(tmax,breaks = 0:horizon, probability = TRUE, 
            main=paste("Timing of peak incidence\nbucket size=",bucket.size),
            xlab="Bucket number", col="purple")
  px = h$density*diff(h$breaks)
  proba.tmax<- data.frame(week=0:horizon, proba=c(px,NA))
  write.csv(proba.tmax,file=paste0(folder,"peakTime_",filesuffix,".csv"))
  
  ### Generation interval (backward GI)
  mean.GI <- GIbck.summary(sim)[["GIbck_mean"]]
  write.csv(mean.GI,file = paste0(folder,"GIbck_",filesuffix,".csv"))
  
}
