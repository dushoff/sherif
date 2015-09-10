library(dplyr)
library(ggplot2)
library(gridExtra)

### Read all output files with max likelihood
###
flist <- system("ls calib_results_maxLik_*.csv", intern = TRUE)
n<- length(flist)

### Merge everything in one dataframe
ff <- list()
for(i in 1:n){
  ff[[i]] <- read.csv(flist[i],header = F)
  names(ff[[i]])<-c("name","value")
  ff[[i]]$statType <- i
}
df <- dplyr::rbind_all(ff)
df2 <- df[(as.character(df$name)!="maxsynlik.synlik"),]

### Distance from true (simulated) values
###
dist <- vector(length=n)
for(i in 1:n){
  dist[i] <- sum(df2$value[df2$statType==i]^2)
}
df.dist <- data.frame(statType=c(1:n), dist=dist)

### PLOTS
###
g.dist <- ggplot(df.dist)+geom_bar(aes(x=factor(statType),
                                       y=dist,
                                       fill=factor(statType)),
                                   stat="identity")

g.prm <- ggplot(df2)+geom_bar(aes(x=name,y=value,fill=factor(statType)),
                              width=0.6,
                              position="dodge",
                              stat="identity")
g.prm <- g.prm+coord_flip()+ylab("Relative Error")+xlab("")

pdf("compstats.pdf",width = 14,height = 7)
grid.arrange(ncol=2, g.dist,g.prm)
dev.off()
