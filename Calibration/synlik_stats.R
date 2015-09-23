########################################################################
###
###   DEFINITION OF SUMMARY STATS FOR SYNTHETIC LIKELIHOOD
###
###   Created 2015-08-11 By David Champredon
###
###   Edited 2015-09-14 : add spatial stats (DC)
###
########################################################################



### Utils functions

norm.reg.coef <- function(x,y,deg){
  # Euclidian norm of polynomial regression coefficients vector
  poly.reg <- lm(formula = y ~ poly(x,deg))
  return(sum(poly.reg$coefficients^2))
}

### ================



### Statistics on the model outputs
### to inform the synthetic likelihood
###
sherif_stats_1 <- function(x,extraArgs,...){
  
  # latest time of the observed data:
  horizon <- extraArgs[["horizon"]]
  
  ts <- x[["ts"]]
  GIbck <- x[["GIbck"]]
  nsim = length(ts)
  
  # All stats
  s1 <- vector(length = nsim)
  s2 <- vector(length = nsim)
  s3 <- vector(length = nsim)
  s4 <- vector(length = nsim)
  s5 <- vector(length = nsim)
  
  for(i in 1:nsim){
    
    # look only up to the largest time 
    # of the observed data (do not take into account of what's beyond)
    df <- subset(ts[[i]], time<horizon)
    
    tt <- df$time
    cumi <- df$cumInc
    cumi_hcw <- df$cumInc_hcw
    cumd <- df$deathsCum
    
    df.inc <- data.frame(t=tt, cuminc=cumi)
    df.death <- data.frame(t=tt, cumd=cumd)
    
    # Cum. incidence general population:
    fit.inc <- glm(formula = cuminc~t, family = "poisson", data= df.inc)
    s1[i] <- fit.inc$coefficients[1]
    s2[i] <- fit.inc$coefficients[2]
    
    # Cum. incidence HCW:
    fit.inc <- glm(formula = cumi_hcw~t, family = "poisson", data= df.inc)
    s6[i] <- fit.inc$coefficients[1]
    s7[i] <- fit.inc$coefficients[2]	
    
    # Cum. deaths:
    fit.death <- glm(formula = cumd~t, family = "poisson", data= df.death)
    s3[i] <- fit.death$coefficients[1]
    s4[i] <- fit.death$coefficients[2]
    
    # Bckwd generation interval
    s5[i] <- mean(GIbck[[i]])
    
    # DEBUG
    # print(s1[i])
    # print(s2[i])
  }
  return(cbind(s1,s2,s3,s4,s5,s6,s7)) #,s8,s9))
}



sherif_stats_2 <- function(x,extraArgs,...){
  
  # latest time of the observed data:
  horizon <- extraArgs[["horizon"]]
  
  ts <- x[["ts"]]
  GIbck <- x[["GIbck"]]
  nsim = length(ts)
  
  # All stats
  s1 <- vector(length = nsim)
  s2 <- vector(length = nsim)
  s3 <- vector(length = nsim)
  s4 <- vector(length = nsim)
  s5 <- vector(length = nsim)
  s6 <- vector(length = nsim)
  
  for(i in 1:nsim){
    
    # look only up to the largest time 
    # of the observed data (do not take into account of what's beyond)
    df <- subset(ts[[i]], time<horizon)
    
    tt <- df$time
    cumi <- df$cumInc
    cumi_hcw <- df$cumInc_hcw
    cumd <- df$deathsCum
    
    # Cum. incidence
    s1[i] <- mean(cumi)	
    s2[i] <- norm.reg.coef(x=tt, y=cumi, deg=3)
    
    # Cum. deaths
    #     s3[i] <- mean(cumd)
    #     s4[i] <- norm.reg.coef(x=tt, y=cumd, deg=3)
    #     
    # Bckwd generation interval
    s3[i] <- mean(GIbck[[i]])
  }
  return(cbind(s1,s2,s3))#,s4,s5))
}



sherif_stats_3 <- function(x,extraArgs,...){
  
  ts <- x[["ts"]]
  GIbck <- x[["GIbck"]]
  nsim = length(ts)
  
  # All stats
  s1 <- vector(length = nsim)
  s2 <- vector(length = nsim)
  s3 <- vector(length = nsim)
  s4 <- vector(length = nsim)
  s5 <- vector(length = nsim)
  s6 <- vector(length = nsim)
  
  for(i in 1:nsim){
    tt <- ts[[i]]$time
    cumi <- ts[[i]]$cumInc
    cumi_hcw <- ts[[i]]$cumInc_hcw
    cumd <- ts[[i]]$deathsCum
    
    # Cum. incidence
    s1[i] <- mean(cumi)	
    
  }
  return(cbind(s1))
}



select.synlik.stats <- function(x){
  
  if(x==1) return(sherif_stats_1)
  if(x==2) return(sherif_stats_2)
  if(x==3) return(sherif_stats_3)
  if(x==4) return(sherif_stats_4)
}

########################################################################
########################################################################
###
###   STATS FOR SPATIAL MODEL
###
########################################################################
########################################################################



sherif_spatial_stats_1 <- function(x,extraArgs,...){
	
	# latest time of the observed data:
	horizon <- extraArgs[["horizon"]]
	
	nLocations <- length(x)
	
	for(loc in 1:nLocations){
		
		ts <- x[[loc]][["ts"]]
		time.firstCase <- x[[loc]][["time_firstCase"]]
		
		nsim = length(ts)
		
		# All stats
		s1 <- vector(length = nsim)
		s2 <- vector(length = nsim)
		s3 <- vector(length = nsim)
		s4 <- vector(length = nsim)
		s5 <- vector(length = nsim)
		s6 <- vector(length = nsim)
		s.gi <- vector(length = nsim)
		
		for(i in 1:nsim){
			# look only up to the largest time 
			# of the observed data (do not take into account of what's beyond)
			df <- subset(ts[[i]], time<horizon)
			
			tt <- df$time
			cumi <- df$cumInc
			cumi_hcw <- df$cumInc_hcw
			cumd <- df$deathsCum
			
			df.inc <- data.frame(t=tt, cuminc=cumi)
			df.death <- data.frame(t=tt, cumd=cumd)
			
			# Cum. incidence general population:
			fit.inc <- glm(formula = cuminc~t, family = "poisson", data= df.inc)
			s1[i] <- fit.inc$coefficients[1]
			s2[i] <- fit.inc$coefficients[2]
			
			if(F){
			  pp <- predict(fit.inc,data.frame(t=tt),type="response")
			  plot(tt,cumi,typ="s",lwd=3, main = paste(round(s1[i],3),round(s2[i],3)))
			  lines(tt,pp,col="red", lwd=3)
			}
			
# 			# Cum. deaths:
# 			fit.death <- glm(formula = cumd~t, family = "poisson", data= df.death)
# 			s3[i] <- fit.death$coefficients[1]
# 			s4[i] <- fit.death$coefficients[2]
			
			# Time first case
			#s6[i] <- mean(time.firstCase[[i]])
		}
		
		Mtmp <- cbind(s1,s2) #,s3,s4)  #,s6)
		if(loc==1) M <- Mtmp
		if(loc>1) M <- cbind(M,Mtmp)
		
		nc <- ncol(M)
		nstats <- ncol(Mtmp)
		tt <- paste0("s",c(1:nstats))
		colnames(M)[(nc-nstats+1):nc] <- paste0("loc",loc,"_",tt)
	}
	
	# Bckwd generation interval
	#
	
	for(i in 1:nsim){
	  m <- numeric(nLocations)
	  for(loc in 1:nLocations){
	    # the GI data are average across all locations
	    GIbck <- x[[loc]][["GIbck"]]
	    m[loc] <- mean(GIbck[[i]])
	  }
	  s.gi[i] <- mean(m,na.rm = T) # <-- remove NAs is important (some epidemic may not have started in some locations, so no GI data!)
# 	  print("==DEBUG::")
# 	  print(m)
# 	  print(s.gi[i])
	}
	
	M <- cbind(M,s.gi)
	colnames(M)[ncol(M)]<- "GI"
	
	return(M)
}



sherif_spatial_stats_2 <- function(x,extraArgs,...){
  
	# latest time of the observed data:
	horizon <- extraArgs[["horizon"]]
	
	nLocations <- length(x)
	
	for(loc in 1:nLocations){
		
		ts <- x[[loc]][["ts"]]
		time.firstCase <- x[[loc]][["time_firstCase"]]
		
		nsim = length(ts)
		
		# All stats
		s1 <- vector(length = nsim)
		s2 <- vector(length = nsim)
		s3 <- vector(length = nsim)
		s4 <- vector(length = nsim)
		s5 <- vector(length = nsim)
		s6 <- vector(length = nsim)
		s.gi <- vector(length = nsim)
		
		for(i in 1:nsim){
			# look only up to the largest time 
			# of the observed data (do not take into account of what's beyond)
			df <- subset(ts[[i]], time<horizon)
			
			tt <- df$time
			cumi <- df$cumInc
			cumi_hcw <- df$cumInc_hcw
			cumd <- df$deathsCum
			
			df.inc <- data.frame(t=tt, cuminc=cumi)
			df.death <- data.frame(t=tt, cumd=cumd)
			
			# Cum. incidence general population:
# 			fit.inc <- glm(formula = cuminc~t, family = "poisson", data= df.inc)
# 			s1[i] <- fit.inc$coefficients[1]
# 			s2[i] <- fit.inc$coefficients[2]
			s1[i] <- mean(cumi)
			s2[i] <- cumi[length(cumi)]
			
		
			
			# 			# Cum. deaths:
			# 			fit.death <- glm(formula = cumd~t, family = "poisson", data= df.death)
			# 			s3[i] <- fit.death$coefficients[1]
			# 			s4[i] <- fit.death$coefficients[2]
			
			# Time first case
			#s6[i] <- mean(time.firstCase[[i]])
		}
		
		Mtmp <- cbind(s1,s2) #,s3,s4)  #,s6)
		if(loc==1) M <- Mtmp
		if(loc>1) M <- cbind(M,Mtmp)
		
		nc <- ncol(M)
		nstats <- ncol(Mtmp)
		tt <- paste0("s",c(1:nstats))
		colnames(M)[(nc-nstats+1):nc] <- paste0("loc",loc,"_",tt)
	}
	
	# Bckwd generation interval
	#
	
	for(i in 1:nsim){
		m <- numeric(nLocations)
		for(loc in 1:nLocations){
			# the GI data are average across all locations
			GIbck <- x[[loc]][["GIbck"]]
			m[loc] <- mean(GIbck[[i]])
		}
		s.gi[i] <- mean(m,na.rm = T) # <-- remove NAs is important (some epidemic may not have started in some locations, so no GI data!)
		# 	  print("==DEBUG::")
		# 	  print(m)
		# 	  print(s.gi[i])
	}
	
	M <- cbind(M,s.gi)
	colnames(M)[ncol(M)]<- "GI"
	
	return(M)
}



select.synlik.spatial.stats <- function(x){
	
	if(x==1) return(sherif_spatial_stats_1)
	if(x==2) return(sherif_spatial_stats_2)
	if(x==3) return(sherif_spatial_stats_3)
	if(x==4) return(sherif_spatial_stats_4)
}

select.spatial.stats <- function(x){
	
	if(x==1) return(sherif_spatial_stats_1)
	if(x==2) return(sherif_spatial_stats_2)
	if(x==3) return(sherif_spatial_stats_3)
	if(x==4) return(sherif_spatial_stats_4)
}