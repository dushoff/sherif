library(plyr)
library(reshape2)
library(snowfall)


run.sherif.parallel <- function(prm.simul, prm.model, ncpus, path.sherif.lib){
	
	sfInit(parallel = FALSE, cpu = ncpus)
	sfLibrary(sherif,lib.loc = path.sherif.lib)
	
	snow.wrap.sherif <- function(i){
		prm.simul.i <- prm.simul
		
		# make sure the seed is different:
		prm.simul.i[["seed"]] <- prm.simul[["seed"]]+i
		# split the number of MC iterations across all cpus:
		prm.simul.i[["mc_iter"]] <- ceiling(prm.simul[["mc_iter"]]/ncpus)
		
		return(rcpp_sherif(paramsSimul = prm.simul.i, 
						   paramsModel = prm.model))
	}
	
	sfExportAll()
	res <- sfSapply(1:ncpus, snow.wrap.sherif,simplify = F)
	sfStop()
	
	n <- length(res)
	res2 <- list()
	for(i in 1:n) {
		if(i==1) res2<- res[[i]]
		if(i>1) res2 = mapply(c, res2,res[[i]],SIMPLIFY = F)
	}
	return(res2)
}


plot.sim.all.mc <- function(sim, varname,...){
	### Plot all MC iterations
	plot(x=sim$time[[1]], y=sim[[varname]][[1]],
		 ylim=c(0,max(unlist(sim[[varname]]))),
		 main = varname, xlab="time",ylab="",
		 typ="s",
		 ...)
	for(i in 2:length(sim$time)) lines(x=sim$time[[i]],y=sim[[varname]][[i]],typ="s",...)
}

agregate.time <- function(sim, varname, timebucket){
	### Aggregate a time-series in coarser time bucket
	###	
	x <- sim[[varname]]
	M <- t(do.call(rbind,x))
	M <- cbind(M,round(sim$time[[1]]/timebucket)*timebucket)
	mc_names <- paste0("mc_",c(1:(ncol(M)-1)))
	colnames(M) <- c(mc_names, "timeround")
	M2 <- melt(as.data.frame(M),measure.vars=mc_names)
	M3 <- ddply(M2,c("variable","timeround"),summarize,n=sum(value))
	tt <- unique(M3$timeround)
	M4 <- as.matrix(dcast(M3,formula = variable ~ timeround,value.var = "n"))
	M5 <- as.matrix(M4[,-1])
	class(M5) <- "numeric"
	return(list(M=M5, time=tt))
}


summarize.sim <- function(sim,varname, timebucket=NULL){
	### Summary of a SHERIF simulation
	###
	if(!is.numeric(timebucket)){
		x <- sim[[varname]]
		M <- do.call(rbind,x)
		tt = sim$time[[1]] ### <-- TO DO: ADD WARNING BC ASSUME ALL TIME VECTORS ARE THE SAME ACROSS MC
	}
	if(is.numeric(timebucket)) {
		a <- agregate.time(sim, varname, timebucket)
		M <- a[["M"]]
		tt <- a[["time"]]
	}
	M.mean <- apply(M,2,FUN=mean)
	M.q2.5 <- apply(M,2,FUN=quantile, probs=0.025)
	M.q10 <- apply(M,2,FUN=quantile, probs=0.10)
	M.q50 <- apply(M,2,FUN=quantile, probs=0.50)
	M.q90 <- apply(M,2,FUN=quantile, probs=0.90)
	M.q97.5 <- apply(M,2,FUN=quantile, probs=0.975)
	
	return(list(time = tt,
				mean = M.mean,
				q2.5 = M.q2.5,
				q10 = M.q10,
				q50 = M.q50,
				q90 = M.q90,
				q97.5 = M.q97.5)
	)
}



GIbck.summary <- function(sim, do.plot = FALSE){
	### Summarize the backward generation interval
	### across all Monte Carlo iterations
	###
	
	# calendar time in the simulation
	# where the nckwd GI was calculated
	t.gi <- sim$GIbck_time[[1]]
	
	# all generaion intervals
	g0 = sim$GIbck_gi
	g = unlist(g0)
	
	if(do.plot) hist(g,breaks = 0:(max(g)+2), 
					 col="lightgrey",
					 main=paste("Histogram of Backward GI at time",round(t.gi)))
	return(list(GIbck_mean=	mean(g),
				GIbck_median=median(g)))
}


plot.summarize.sim <- function(sim, varname, timebucket=NULL, typ, ...){
	### Plot the summary of a SHERIF simulation
	###
	z <- summarize.sim(sim,varname, timebucket)
	ymax <- max(z[["q97.5"]])
	tt <- z[["time"]]
	
	plot(x=tt, y=z[["mean"]],
		 lwd = 6, xlab="time", ylab=varname, las=1,
		 ylim=c(0,ymax), typ=typ,
		 ...)
	lines(x=tt, y=z[["q2.5"]], typ=typ,lwd=1, col="lightgrey")
	lines(x=tt, y=z[["q97.5"]], typ=typ,lwd=1, col="lightgrey")
	
	lines(x=tt, y=z[["q10"]], typ=typ, lwd=2, col="grey")
	lines(x=tt, y=z[["q90"]], typ=typ,lwd=2, col="grey")
}


plot.Reff <- function(sim, timebucket){
	### Plot Effective reproductive number
	### at different dates and 
	### averaged across all simulations
	
	n <- length(sim$Reff_timeAcq)
	
	for(i in 1:n){
		tacq <- sim$Reff_timeAcq[[i]]
		tacq2 <- round(tacq/timebucket)*timebucket
		nseck <- sim$Reff_n2ndCases[[i]]
		df <- data.frame(tacq2=tacq2, nseck=nseck, mc=i)
		df2 <- ddply(df,c("tacq2","mc"),summarize,Reff=mean(nseck)) 
		if(i==1) X = df2
		if(i>1) X = rbind(X,df2)
	}
	X2 <- ddply(X,"tacq2",summarize,m=mean(Reff),min=min(Reff),max=max(Reff))
	library(ggplot2)
	g <- ggplot(X2) + geom_linerange(aes(x=tacq2,y=m,ymin=min,ymax=max)) 
	g <- g + geom_point(aes(x=tacq2,y=m))
	g <- g + xlab("time")+ylab("Reff (avg)") + ggtitle("Effective Reproductive Number")
	plot(g)
}

