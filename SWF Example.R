

source(paste(getwd(),"SWF Algorithm.R",sep="/"))
library(ggridges)
library(ggplot2)

# initial parameters
my_theta <- c(0.5,0.5)
# data -> from binomial n=10
y <- c(1,1,0,1)

# update pi0
x1 <- seq(0.01,1-0.01,by=0.01)
newpars <- c(my_theta[1]+y[1],my_theta[2]+1-y[1])
plot(dbeta(x=x1,shape1=newpars[1],shape2=newpars[2]),type="l",col="black",lwd=2)
lines(dbeta(x=x1,shape1=my_theta[1],shape2=my_theta[2]),type="l",col="black",lwd=2)
# predict pi1
newpars
out <- swf.weights(y[2],1,2,newpars)
# update pi1
newpars1 <- c(NA,NA)
newpars1[1] <- newpars[1]+y[2]
newpars1[2] <- newpars[2]+1-y[2]
pi1_dens <- out[1,1]*dbeta(x1,newpars1[1]+0,newpars1[2])+out[2,1]*dbeta(x1,newpars1[1]+1,newpars1[2])
lines(pi1_dens,type="l",lwd=2,col="purple")
# predict pi2
newpars1
out <- swf.weights(y[3],1,2,newpars1)
# update pi2
newpars1 <- c(NA,NA)
newpars1[1] <- newpars[1]+y[3]
newpars1[2] <- newpars[2]+1-y[3]
pi2_dens <- pi1_dens*( out[1,1]*dbeta(x1,newpars1[1],newpars1[2]+0)+out[1,2]*dbeta(x1,newpars1[1],newpars1[2]+1) ) 
lines(pi2_dens,type="l",lwd=2,col="pink")





lines(out[[2]],col="green",type="l",lwd=2)
out2 <- swf.weights(0,1,2,c(0.5,0.5))
lines(out2[[2]],col="red",type="l",lwd=2)
theta <- c(0.5,0.5)

plot(0.86*dbeta(x1,theta[1]+0,theta[2])+0.14*dbeta(x1,theta[1]+1,theta[2]),type="l",lwd=2,col="green")
lines(0.86*dbeta(x1,theta[1],theta[2]+0)+0.14*dbeta(x1,theta[1],theta[2]+1),type="l",lwd=2,col="red")

# not quite mixing properly


# surely we should have a different distribution here?


x1 <- seq(0.01,1-0.01,by=0.01)
plot(density(dbeta(x=x1,shape1=0.5,shape2=0.5)),ylab="Density",ylim=c(0,10),type="l",col="blue")
out <- swf.weights(1,1,2,c(0.5,0.5))
lines(out[[2]],col="orange")
out2 <- swf.weights(0,1,2,c(0.5,0.5))
lines(out2[[2]],col="red")
#plot(dbeta(x=x1,shape1=0.5+3+3,shape2=0.5+6+2),ylab="Density",ylim=c(0,10),type="l",col="blue")

out2 <- swf.weights(17,20,100,c(0.5,0.5))
plot(out2[[2]])



