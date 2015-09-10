###
###    ANALYZE RESULTS FROM RUNNING CALIBRATION
###    SEVERAL TIMES WITH SAME TARGET PARAMATERS
###    BUT DIFFERENT RANDOM NUMBER GENERATOR SEEDS
###

library(ggplot2)
library(dplyr)

# Retrieve true parameters:
prm.true <- read.csv("param_true4fit.csv",header = F)

# Name of the parameter that were fitted:
prm.name <- c("beta_IS","meanDur_latent")

# Retrieve results:
fls <- system("ls calib_results_post_?_*.csv", intern = T)
n <- length(fls)
flist <- list()
for(i in 1:n) {
  flist[[i]] <- read.csv(fls[i],header = T)
  names(flist[[i]])[1] <- c("param")
  flist[[i]]$seed <- i
}
df <- dplyr::rbind_all(flist)

# add the true value information to the results data frame:
df$trueval <- NA
for(i in 1:length(prm.name)) 
  df$trueval[df$param==prm.name[i]]<- prm.true[prm.true$V1==prm.name[i],2]

all.res <- subset(df, param %in% prm.name)

pdf("analyze_multiseed.pdf",width=15,height = 10)

g <- ggplot(all.res)+geom_boxplot(aes(x=param,y=post.mean))+facet_wrap(~param,scales = "free_y")
g <- g + geom_point(aes(x=param,y=trueval),colour="red",size=3)
g <- g + ggtitle("posterior mean")
print(g)

g <- ggplot(all.res)+geom_pointrange(aes(x=seed,y=post.mean,ymin=q2.5, ymax=q97.5))
g <- g + geom_pointrange(aes(x=seed,y=post.mean,ymin=q10, ymax=q90),size=1.5)
g <- g + geom_hline(aes(yintercept=trueval),colour="red",size=2) 
g <- g + facet_wrap(~param,scales = "free_y")
g <- g + ggtitle("posterior mean ; 80 & 95% CI")
plot(g)


dev.off()