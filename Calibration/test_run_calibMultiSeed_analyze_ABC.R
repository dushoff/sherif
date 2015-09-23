
system("rm best_ABC*.csv",intern = F, wait = T)

n <- 10
for(i in 1:n){
  cmd <- paste0("Rscript calibration_spatial_ABC_genData.R 1 ",i)
  system(command = cmd, intern = F, wait = TRUE)
}

# Retrieve true parameters:
prm.true <- as.data.frame(read.csv("true_prm.csv",header = T)); prm.true

fls <- system("ls best_ABC_*.csv", intern = T, wait = T)
flist <- list()
for(i in 1:n) {
  flist[[i]] <- read.csv(fls[i],header = T)
  # names(flist[[i]])[1] <- c("param")
  flist[[i]]$seed <- i
}
df <- dplyr::rbind_all(flist) ; df


summary(df)
library(reshape2)
library(ggplot2)

df.m <- melt(data = df)
df.m <- subset(df.m, variable!="ABCdist")
df.m <- subset(df.m, variable!="seed")

df.m$trueval <- NA
nn <- length(names(prm.true))
for(k in 1:nn){
  tmp <- names(prm.true)[k]
  df.m$trueval[df.m$variable==tmp] <- prm.true[1,tmp]
}

pdf("multi_seed_test.pdf", width=15,height = 10)
g <- ggplot(df.m)+geom_boxplot(aes(x=variable,y=value),alpha=0.4) +geom_point(aes(x=variable,y=value)) 
g <- g + geom_point(aes(x=variable,y=trueval),size=6,colour="red")
g <- g + geom_violin(aes(x=variable,y=value),alpha=0.2) 
g <- g + facet_wrap(~variable, scales="free_y")+ggtitle("Multiseed Test")+xlab("")+ylab("")
plot(g)
dev.off()

# Name of the parameter that were fitted:
prm.name <- c("beta_IS","meanDur_latent")
