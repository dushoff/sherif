norm.reg.coef <- function(x,y,deg){
  # Euclidian norm of polynomial regression coefficients vector
  poly.reg <- lm(formula = y ~ poly(x,deg))
  return(sum(poly.reg$coefficients^2))
}


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


#############################################
#### OLD STUFF TO KEEP FOR REFERENCE....

# sherif_stats_2 <- function(x,extraArgs,...){
#   
#   ts <- x[["ts"]]
#   GIbck <- x[["GIbck"]]
#   nsim = length(ts)
#   
#   # All stats
#   s1 <- vector(length = nsim)
#   s2 <- vector(length = nsim)
#   s3 <- vector(length = nsim)
#   s4 <- vector(length = nsim)
#   s5 <- vector(length = nsim)
#   s6 <- vector(length = nsim)
#   s7 <- vector(length = nsim)
#   s8 <- vector(length = nsim)
#   s9 <- vector(length = nsim)
#   
#   for(i in 1:nsim){
#     tt <- ts[[i]]$time
#     cumi <- ts[[i]]$cumInc
#     cumi_hcw <- ts[[i]]$cumInc_hcw
#     cumd <- ts[[i]]$deathsCum
#     
#     # Cum. incidence
#     x <- lm(cumi~poly(tt,degree = 3))$coefficients
#     s1[i] <- x[1] #mean(cumi)	
#     s2[i] <- x[2] #norm.reg.coef(x=tt, y=cumi, deg=3)
#     s3[i] <- x[3]
#     s4[i] <- x[4]
#     
#     # Cum. incidence HCW
#     # s6[i] <- mean(cumi_hcw)	
#     
#     # Cum. deaths
#     # s3[i] <- norm.reg.coef(x=tt, y=cumd, deg=3)
#     x <- lm(cumd~poly(tt,degree = 3))$coefficients
#     s5[i] <- x[1] #mean(cumi)	
#     s6[i] <- x[2] #norm.reg.coef(x=tt, y=cumi, deg=3)
#     s7[i] <- x[3]
#     s8[i] <- x[4]
#     
#     # Bckwd generation interval
#     s9[i] <- mean(GIbck[[i]])
#     # s5[i] <- sd(GIbck[[i]])
#   }
#   return(cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9))
# }
