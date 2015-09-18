### R package for SHERIF is installed locally:
path.sherif.lib <- "../lib"

library(sherif,lib.loc = path.sherif.lib)
library(synlik)

source("loadParam_FCT.R")
source("wrap_FCT.R")

### Load the 'true' parameters that will
### generate the target data 
###
prm.model.true <- loadParam("param_true4fit.csv")
prm.simul.true <- loadParam("param_simul_true4fit.csv")

### Load parameters to thin simulated data
prm.gen.data <- read.csv("generate-target-data-param.csv")

n <- 5
par(mfrow=c(n,n))
cn <- vector(length = n^2)
cn2 <- vector(length = n^2)
cn3 <- vector(length = n^2)
cn4 <- vector(length = n^2)

for(i in 1:n^2){
  set.seed(i*13)
  ### simulate one iteration to generate target data
  ###
  true.data <- sherif_wrap(param = c(beta_IS = prm.model.true[["beta_IS"]]), # <-- does not matter which param chosen here. Just to comply with function required signature
                           nsim = 1, 
                           extraArgs = list(
                             "paramsSimul" = prm.simul.true, 
                             "paramsModel" = prm.model.true))
  
  
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
  
  ###############################################################################################
  ###############################################################################################
  
  ts.plot <- function(ts,var, title,...){
    head(ts)
    typ <- "s"
    if(length(ts$time)>50) typ <- "l"
    plot(x=ts$time, y=ts[,var], typ=typ, lwd=6, 
         main=title, xlab="time", ylab=var,...)
    grid()
  }
  
  ts.true <- true.data2[["ts"]][[1]]
  ts.plot(ts.true,"cumInc", "Target Data - Cum Inc",ylim=c(0,500))
  
  tt <- ts.true$time
  cumi <- ts.true$cumInc
  
  lm1 <- lm(cumi~poly(tt,degree = 3))
  
  zz <- data.frame(tt=seq(0,max(tt),length.out = 100))
  yy <- predict(lm1,zz)
  lines(x=zz$tt, y=yy, lwd=3,col="red")
  cn[i] <- sum(lm1$coefficients^2)
  cn2[i]<-cumi[length(cumi)]
  cn3[i]<-cumi[length(cumi)]-cumi[round(length(cumi)/2)]
  cn4[i]<-mean(cumi)
  
}

sd(cn)/mean(cn)
sd(cn2)/mean(cn2)
sd(cn3)/mean(cn3)
sd(cn4)/mean(cn4)
