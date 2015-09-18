plot.pair.synlik  <- function(df,var1, var2, max.synlik){
  
  if(var1!="synlik" & var2!="synlik"){
    lk <- df$synlik
    rlk <- rank(lk)
    df$cex <- rlk/length(lk)*2
    col <- rgb(0,0,1,alpha = 0.7)
    
    # Results
    plot(x=df[,var1],y=df[,var2],
         pch=15, cex= df$cex ,
         col=col,
         xlab=var1, ylab=var2)
    points(x=max.synlik[,var1], y=max.synlik[,var2], 
           col="red", cex=3, lwd=3)
    
    # True values (if provided)
    if(!is.null(true.prm)){
      var1.true <- as.numeric(true.prm[var1])
      var2.true <- as.numeric(true.prm[var2])
      print(var1.true); print(var2.true)
      abline(v=var1.true, lty=2, col="orange",lwd=2)
      abline(h=var2.true, lty=2, col="orange",lwd=2)
    }
  }
}



inference_LHS <- function(prm2fit.val,
                          n.cpu,
                          nsim,
                          path.sherif.lib,
                          do.plot,
                          true.prm=NULL){
  
  ### LATIN HYPERCUBE SAMPLING 
  ### TO FIND MAXIMUM SYNTHETIC LIKELIHOOD
  
  # Integrity checks:
  stopifnot(length(prm2fit.val)>0)
  stopifnot(nsim>4)
  stopifnot(n.cpu>0)
  if(!is.null(true.prm)) stopifnot(names(prm2fit.val)==names(true.prm))
  
  # Snowfall setup for parallel computing:
  sfInit(parallel = TRUE, cpu = n.cpu)
  sfLibrary(sherif, lib.loc = path.sherif.lib)
  sfLibrary(synlik)
  
  sherif_snowWrap <- function(i){
    ###
    ### WRAP FOR USE WITH 'SNOWFALL' PACKAGE
    ###
    param.i <- as.numeric(prm2fit.val[i,])
    names(param.i) <- names(prm2fit.val)
    message(paste("debug:: ",param.i)) # <-- DEBUG
    res <- slik(object = sherif_sl, 
                param = param.i, 
                multicore = F,  nsim=nsim)
    return(res)
  }
  
  idx.apply <- 1:npts.lhs
  
  ### Parallel execution:
  sfExportAll()
  system.time(xx <- sfSapply(idx.apply, sherif_snowWrap))
  sfStop()
  
  ### Data frame holding all values computed:
  df <- cbind(prm2fit.val, synlik = xx)
  
  ### Maximum synthetic likelihood:
  max.synlik <- df[which.max(df$synlik),]
  
  if(do.plot){
    np <- length(names(df))
    pp <- floor(sqrt(np*(np+1)/2))
    par(mfrow=c(pp,pp))
    for(i in 1:(np-2)) 
      for(j in (i+1):(np-1)) 
        plot.pair.synlik(df,names(df)[i], names(df)[j], max.synlik)
    for(i in 1:(np-1)) 
      if(names(df)[i]!="synlik")
        plot(x=df[,i], y=exp(df$synlik),
             pch=16,xlab=names(df)[i],ylab="EXP(synlik)")
  }
  return(list(bestfit = max.synlik, alltrials = df ))
}



inference_spatial_LHS <- function(prm2fit.val,
						  n.cpu,
						  nsim,
						  path.sherif.lib,
						  do.plot,
						  true.prm=NULL){
	
	### LATIN HYPERCUBE SAMPLING 
	### TO FIND MAXIMUM SYNTHETIC LIKELIHOOD
	### FOR SPATIAL MODEL
	
	# Integrity checks:
	stopifnot(length(prm2fit.val)>0)
	stopifnot(nsim>4)
	stopifnot(n.cpu>0)
	#if(!is.null(true.prm)) stopifnot(names(prm2fit.val)==names(true.prm))
	
	# Snowfall setup for parallel computing:
	sfInit(parallel = TRUE, cpu = n.cpu)
	sfLibrary(sherif, lib.loc = path.sherif.lib)
	sfLibrary(synlik)
	
	sherif_snowWrap <- function(i){

		### WRAP FOR USE WITH 'SNOWFALL' PACKAGE

		param.i <- as.numeric(prm2fit.val[i,])
		names(param.i) <- names(prm2fit.val)
		
		print(param.i) # <-- DEBUG
		
		res <- try(slik(object = sherif_spatial_sl, 
		                param = param.i, 
		                multicore = F,  nsim=nsim), 
		           silent = TRUE)
		
		# * * * IMPORTANT NOTE * * * 
		#
		# if the synthetic likelihood fails, it is set to 0 (i.e. log lik large negative)
		# (failure is likely to happen when a summary stat
		# is NA or constant[:empirical variance=0 => divide by 0 => lik=Inf])
		if (class(res)=="try-error") res <- -9E99
		
		return(res)
	}
	
	idx.apply <- 1:npts.lhs
	
	### Parallel execution:
	sfExportAll()
	system.time(xx <- sfSapply(idx.apply, sherif_snowWrap))
	sfStop()
	
	### Data frame holding all values computed:
	df <- cbind(prm2fit.val, synlik = xx)
	
	### Maximum synthetic likelihood:
	max.synlik <- df[which.max(df$synlik),]
	
	if(do.plot){
		np <- length(names(df))
		pp <- floor(sqrt(np*(np+1)/2))
		par(mfrow=c(pp,pp))
		for(i in 1:(np-2)) 
			for(j in (i+1):(np-1)) 
				plot.pair.synlik(df,names(df)[i], names(df)[j], max.synlik)
		for(i in 1:(np-1)) 
			if(names(df)[i]!="synlik")
				plot(x=df[,i], y=exp(df$synlik),
					 pch=16,xlab=names(df)[i],ylab="EXP(synlik)")
	}
	return(list(bestfit = max.synlik, alltrials = df ))
}







inference_mcmc <- function(synlik.object,
                           prm.init,
                           mcmc.iter, mcmc.warmup,
                           nsim,
                           do.plot){
  ###   MCMC inference
  ###
  
  n <- length(prm.init)
  
  mcmc_sherif_sl <- smcmc(synlik.object, 
                          initPar = prm.init,
                          niter = mcmc.iter, 
                          burn = mcmc.warmup,
                          priorFun = function(input, ...) sum(input), 
                          propCov = diag(c(0.01,0.01,0.01)), #diag(x = rep(sum(prm.init^2),n)/100), 
                          nsim = nsim)
  
  res <- list()
  for(i in 1:n) res[[i]] <- mcmc_sherif_sl@chains[,i]
  
  if(do.plot){
    par(mfrow=c(n,2))
    for(i in 1:n) {
      plot(mcmc_sherif_sl@chains[,i])
      hist(mcmc_sherif_sl@chains[,i], breaks=30, col="grey")
    }
  }
  
}



diff.from.true <- function(true.prm,max.synlik, rel = TRUE){
  ### DIFFERENCE B/W CALIBRATION RESULTS AND TRUE PARAMETER VALUES
  ### 
  mx <- unlist(c(max.synlik[names(true.prm)]))
  trueprm <- unlist(true.prm)
  d <- mx-trueprm
  if(rel) d<- (mx-trueprm)/trueprm
  return(unlist(c(d,maxsynlik=max.synlik["synlik"])))
}



posteriors <- function(df, true.prm=NULL, do.plot=TRUE, prior.range=NULL){
  
  ### Keep only the top [25]% likelihood
  #DEBUG
  print(str(df))
  print(df)
  lik.thresh <- quantile(df$synlik,probs=0.75, na.rm = T)
  df <- subset(df, synlik>= lik.thresh)
  
  ### Stats for posterior
  post.mean <- apply(df, MARGIN = 2,FUN = mean, na.rm = T)
  post.median <- apply(df, MARGIN = 2,FUN = median, na.rm = T)
  q2.5 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.025, na.rm = T)
  q10 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.10, na.rm = T)
  q90 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.90, na.rm = T)
  q97.5 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.975, na.rm = T)
  
  ### Store everything in a data frame
  res = data.frame(post.mean)
  res = cbind(res,post.median)
  res = cbind(res,q2.5)
  res = cbind(res,q10)
  res = cbind(res,q90)
  res = cbind(res,q97.5)
  
  if(do.plot){
    n = length(names(df))
    nr = nrow(df)
    nn = ceiling(sqrt(n))
    par(mfrow=c(nn,nn))
    for(i in 1:n) {
      if(names(df)[i]!="synlik") { # <-- don't care about histogram of likelihood
        # Define breaks such that the range of prior
        # (if specified) is displayed in the plot:
        brks <- min(nr,20)
        if(!is.null(prior.range)) {
          brks <- seq(prior.range[1,names(df)[i]],
                      prior.range[2,names(df)[i]],
                      length.out = 20)
        }
        h<-hist(df[,i],main=names(df)[i], 
                xlab="",yaxt="n",
                freq = FALSE,
                col="lightgrey",border="grey",
                breaks = brks)
        h.ymax<- max(h$density)
        idx <- row.names(res)==names(df)[i]
        points(x=res$post.mean[idx],y=0.1*h.ymax, pch=16,cex=3)
        points(x=res$post.median[idx],y=0.1*h.ymax, pch="|",cex=2)
        segments(x0=res$q2.5[idx],y0=0.1*h.ymax, 
                 x1=res$q97.5[idx],y1=0.1*h.ymax,
                 lwd=2)
        segments(x0=res$q10[idx],y0=0.1*h.ymax, 
                 x1=res$q90[idx],y1=0.1*h.ymax,
                 lwd=6)
        if(!is.null(true.prm)){
          try(abline(v=true.prm[names(df)[i]][[1]], lty=2, lwd=3,col="orange"))
        }
      }
    }
  }# end if do.plot
  
  return(res)
}


posteriors__OLD <- function(df, true.prm=NULL, do.plot=TRUE){
  
  ### Keep only the top [25]% likelihood
  lik.thresh <- quantile(df$synlik,probs=0.75)
  df <- subset(df, synlik>= lik.thresh)
  
  ### Stats for posterior
  post.mean <- apply(df, MARGIN = 2,FUN = mean)
  post.median <- apply(df, MARGIN = 2,FUN = median)
  q2.5 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.025)
  q10 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.10)
  q90 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.90)
  q97.5 <-  apply(df, MARGIN = 2,FUN = quantile, probs=0.975)
  
  ### Store everything in a data frame
  res = data.frame(post.mean)
  res = cbind(res,post.median)
  res = cbind(res,q2.5)
  res = cbind(res,q10)
  res = cbind(res,q90)
  res = cbind(res,q97.5)
  
  if(do.plot){
    n = length(names(df))
    nr = nrow(df)
    nn = ceiling(sqrt(n))
    par(mfrow=c(nn,nn))
    for(i in 1:n) {
      h<-hist(df[,i],main=names(df)[i], xlab="",yaxt="n",
              freq = FALSE,
              col="lightgrey",border=NA,
              breaks = min(nr,20))
      h.ymax<- max(h$density)
      idx <- row.names(res)==names(df)[i]
      points(x=res$post.mean[idx],y=0.1*h.ymax, pch=16,cex=3)
      points(x=res$post.median[idx],y=0.1*h.ymax, pch="|",cex=2)
      segments(x0=res$q2.5[idx],y0=0.1*h.ymax, 
               x1=res$q97.5[idx],y1=0.1*h.ymax,
               lwd=2)
      segments(x0=res$q10[idx],y0=0.1*h.ymax, 
               x1=res$q90[idx],y1=0.1*h.ymax,
               lwd=6)
      if(!is.null(true.prm)){
        try(abline(v=true.prm[names(df)[i]][[1]], lty=2, lwd=3,col="orange"))
      }
    }
  }# end if do.plot
  
  return(res)
}


