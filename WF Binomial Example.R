
source(paste(getwd(),"SWF Algorithm.R",sep="/"))
library(ggridges)
library(ggplot2)

# initial parameters
my_theta <- c(1.5,1.5)
# data -> from binomial n=5
nt <- 5
y <- c(3,4,2,2,5)
tdiff <- rep(0.001,5)

pi_dens <- dbeta(x=x1,shape1=my_theta[1],shape2=my_theta[2])

# update pi0
x1 <- seq(0.01,1-0.01,by=0.01)
newpars <- c(my_theta[1]+y[1],my_theta[2]+nt-y[1])
pi0_dens <- dbeta(x=x1,shape1=newpars[1],shape2=newpars[2])



out1 <- Phi(y[2],nt,tdiff[2],newpars,prior="wf")
pi1_dens <- out1[[1]]
newpars1 <- out1[[2]]


out2 <- Phi(y[3],nt,tdiff[3],newpars1,prior="wf")
pi2_dens <- out2[[1]]
newpars2 <- out2[[2]]


out3 <- Phi(y[4],nt,tdiff[4],newpars2,prior="wf")
pi3_dens <- out3[[1]]
newpars3 <- out3[[2]]


out4 <- Phi(y[5],nt,tdiff[5],newpars3,prior="wf")
pi4_dens <- out4[[1]]
newpars4 <- out4[[2]]


library(tidyr)
library(dplyr)
library(reshape2)

nl <- length(pi1_dens)

posteriors <- matrix(c(pi_dens,pi0_dens,pi1_dens,pi2_dens,pi3_dens,pi4_dens),nrow=6,ncol=nl,byrow=T)
rownames(posteriors) <- c("pi","pi0","pi1","pi2","pi3","pi4")


post_melt <- melt(posteriors)
ggplot(data=post_melt,aes(x=value,y=Var1))+geom_density_ridges2()


