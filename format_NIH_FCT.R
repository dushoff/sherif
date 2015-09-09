distrib.output.sherif <- function(sim, varname, 
								  t, interv.proba, 
								  timebucket = NULL,
								  do.plot=TRUE){
	### Extract distribution of cumulative incidence at time 't'
	### (agregation in time bucket possible, e.g. for incidence (not cum inc!))
	###
	
	if(!is.numeric(timebucket)){
		res <- distrib.output.sherif.noAgreg(sim, 
											 varname, 
											 t, 
											 interv.proba, do.plot)	
	}
	if(is.numeric(timebucket)) {
		a <- agregate.time(sim, varname, timebucket)
		M <- a[["M"]]
		tt <- a[["time"]]
		
		idx.t <- which(abs(tt-t)<timebucket/2)
		t.bucket <- tt[idx.t]
		title <- paste0("Distribution of ",varname,
						" at time ",t.bucket, 
						"\n(time ",t," asked, bucket=",timebucket,")")
		
		h <- hist(M[,idx.t],breaks = interv.proba,
				  plot = do.plot,
				  probability = TRUE, col = "grey",
				  main = title, xlab="", las=1)
		
		px = h$density*diff(h$breaks)
		proba.interval <- data.frame(interval=interv.proba, proba=c(px,NA))
		res <- list(hist = h,
					proba.interval = proba.interval)
	}
	return(res)
}


distrib.output.sherif.noAgreg <- function(sim, varname, t, interv.proba, do.plot=TRUE){
	
	### Extract distribution of cumulative incidence at time 't'
	
	mc.iter <- length(sim$time)
	### min and max across all Monte Carlo
	min.val <- min(unlist(sim[[varname]]))
	max.val <- max(unlist(sim[[varname]]))
	
	x <- numeric(mc.iter)
	for(i in 1:mc.iter){
		simtime <- sim$time[[i]]
		dt <- min(diff(simtime))
		idx.t <- which(abs(simtime-t)<dt/2)
		if(length(idx.t)==0){stop(paste0("time ",t," not found in simulation!"))}
		
		x[i] <- sim[[varname]][[i]][idx.t]  
	}
	
	h <- hist(x, probability = TRUE,
			  main = paste(varname,"at time",t),
			  xlab="",
			  breaks = interv.proba,# seq(min(x)-1,max(x)+1,by=1),
			  col="grey", plot = do.plot)
	
	px = h$density*diff(h$breaks)
	proba.interval <- data.frame(interval=interv.proba, proba=c(px,NA))
	return(list(hist=h,
				proba.interval=proba.interval))
}


time.peakIncidence <- function(sim,timebucket){
  ### Calculate time of peak incidence
  ###
  a <- agregate.time(sim = sim,varname = "incidence",timebucket = 7)
  idx.tmax <- apply(a$M,MARGIN = 1, FUN = which.max)
  tmax <- a$time[idx.tmax]
  return(tmax)  
}


