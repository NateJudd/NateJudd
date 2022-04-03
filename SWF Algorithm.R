

# data = vector of 0s and 1s representing observations of the allelle
# init - inital values of the beta distribution


swf <- function(data,times,init){
  
  niter <- length(data)
  alpha0 <- init[1]
  beta0 <- init[2]
  n0 <- 0
  p0 <- 0
  
  # output list
  post <- list()
  
  # initialise
  x1 <- seq(0,1,by=0.01)
  prior <- dbeta(x=x1,shape1=alpha0,shape2=beta0)
  alpha <- rep(NA,niter)
  beta <- rep(NA,niter)
  mix <- x1
  
  # update
  n <- n0 + data[1]
  p <- p0 + (1-data[1])
  alpha[1] <- alpha0 + n
  beta[1] <- beta0 + p
  post[[1]] <- dbeta(x=x1,shape1=alpha[1],shape2=beta[1])
  
  # iterate over the data
  for (n_data in 2:niter){
    
    weights <- swf.weights(i=n,j=p,t=times[n],pars=init)
    print(weights)
    readline()
    
    ### predict
    # initialise mixture distribution object
    
    # iterate over p=0 to n to construct the mixture
    for (p in 1:n){
      # update too
      alpha[n] <- alpha[1] + p + data[n]
      beta[n] <- beta[2] + (n-p) + (1-data[n])
      mix <- mix + weights[n,p] * dbeta(x=x1,shape1=alpha[n],shape2=beta[n])
    }
    
    # update
    
    # save posterior
    post[[n]] <- mix
    # update parameters
    n0 <- n0 + n
    p0 <- p0 + p
  }
  
}

# i,j versus i,j-i

swf <- function(data,times,init){
  # extract initial inputs
  n_datum <- length(data) # vector input - 0s and 1s
  theta0 <- init # alpha and beta for the initial/stationary distribution
  tdiff <- times[1]
  # set up output data storage
  posteriors <- list()
  # update to extract the first posterior
  i <- data[1]
  j <- (1-data[1])
  #theta0 + c(data[1],1-data[1])
  # iterate over the data
  for (datum in 1:n_datum){
    posteriors[[datum]]
    # predict
    swf.predict(i=i,j=j,t=tdiff,pars=theta0)
    # update
    swf.update()
    # update the time parameter for the next iteration
    tdiff <- times[datum+1] - times[datum]
  }
  
  return(posteriors)
}

swf.predict <- function(i,j,t,pars){
  # initialise density sequence
  x1 <- seq(0,1,by=0.01)
  # extract the new weights for \nu_{i,j-i} input
  weights <- swf.weights(i,j,t,pars)
  new_distribution <- weights * dbeta(x=x1,shape1=i+pars[1],shape2=j+pars[2])
}


wf.weights <- function(n_succ,n,t,pars){
  # extratc initial pars
  i <- n_succ
  j <- n
  # how many iterations?
  nw <- 2*(i+j+1)
  theta <- pars
  
  # iterate over the double index k=0,...,i  and l=0,...,j-i
  # initialise output matrix
  nr <- i+1
  nc <- (j-i+1)
  weights <- matrix(NA,nrow=nr,ncol=nc)
  # initialise distribution/posterior
  x1 <- seq(0.01,1-0.01,by=0.01)
  new_distribution <- rep(0,99)
  for (k in 0:i){
    for (l in 0:(j-i)){
      # easy case
      if(k==0 & l==0){
        weights[i-k+1,j-i-l+1] <- exp(- swf.lambda(n=j,pars=theta) * t)
      }else{
        # 3 parts --> hypergeometric prob func, rates product, sum component
        hyp <- choose(i,k)*choose(j-i,l) / choose(j,k+l)
        # initialise rates_prod
        rates_prod <- 1
        # multiply rates together
        for (num in 0:(k+l-1)){
          rates_prod <- rates_prod * swf.lambda(n=j-num,pars=theta)
        }
        # initialise sum_comp
        sum_comp <- 0
        # add on ratios of weights
        # compute B_t(a_j,...,a_{j-k-l})
        for(num1 in 0:(k+l)){
          # compute denominator rate-difference product
          denom <- 1
          for (num2 in 0:(k+l)){
            if(num2 != num1){
              denom <-  denom * abs( swf.lambda(n=j-num1,pars=theta) - swf.lambda(n=j-num2,pars=theta) )
            }else{
              denom <- denom
            }
          }
          sum_comp <- sum_comp + exp(-swf.lambda(n=j-num1,pars=theta)*t)*(-1)^(num1) / denom
        }
        
        weights[i-k+1,j-i-l+1] <- hyp * rates_prod * (-1)^(k+l) * sum_comp
        
        }
      
      }
  }
  # weights should sum to 1
  
  return(weights)
}



swf.weights <- function(n_succ,n,t,pars){
  # extratc initial pars
  i <- n_succ
  j <- n
  # how many iterations?
  nw <- 2*(i+j+1)
  theta <- pars
  
  # iterate over the double index k=0,...,i  and l=0,...,j-i
  # initialise output matrix
  nr <- i+1
  nc <- (j-i+1)
  weights <- matrix(NA,nrow=nr,ncol=nc)
  # initialise distribution/posterior
  x1 <- seq(0.01,1-0.01,by=0.01)
  new_distribution <- rep(0,99)
  for (k in 0:i){
    for (l in 0:(j-i)){
      # easy case
      if(k==0 & l==0){
        weights[i-k+1,j-i-l+1] <- exp( swf.psi(n=j,pars=theta) * t)
      }else{
        # 3 parts --> hypergeometric prob func, rates product, sum component
        hyp <- choose(i,k)*choose(j-i,l) / choose(j,k+l)
        # initialise rates_prod
        rates_prod <- 1
        # multiply rates together
        for (num in 0:(k+l-1)){
          rates_prod <- rates_prod * swf.lambda(n=j-num,pars=theta)
        }
        # initialise sum_comp
        sum_comp <- 0
        # add on ratios of weights
        # compute B_t(a_j,...,a_{j-k-l})
        for(num1 in 0:(k+l)){
          # compute denominator rate-difference product
          denom <- 1
          for (num2 in 0:(k+l)){
            if(num2 != num1){
              denom <-  denom * abs( swf.lambda(n=j-num1,pars=theta) - swf.lambda(n=j-num2,pars=theta) )
            }else{
              denom <- denom
            }
          }
          sum_comp <- sum_comp + exp( swf.psi(n=j-num1,pars=theta) * t ) *(-1)^(num1) / denom
        }
        
        weights[i-k+1,j-i-l+1] <- hyp * rates_prod * (-1)^(k+l) * sum_comp
        
      }
      
    }
  }
  # weights should sum to 1
  
  return(weights)
}


swf.psi <- function(n,pars,rho=1){
  # hard code integration limits
  my_lim <- 4
  # compound poisson process - b_0=0, N with intensity rho, xi N(0,1)-valued RV
  rate <- swf.lambda(n,pars)
  expectation <-  integrate(f=my_integrand,u=rate,lower=-my_lim,upper=my_lim)
  new_par <- rho * expectation$value
  return(new_par)
}

# compound poisson process with b_0=0
my_integrand <- function(z,u=1){
  #z <- dnorm(seq(0, 1, by = 0.01), 0, 1)
  zdens <- dnorm(z,0,1)
  fz <- ( 1-exp(u*zdens) )
  return(fz)
}


swf.lambda <- function(n,pars){
  # WF's dual's rates --> Kingman's Coalescent rates
  rate <- n * ( 2*(n-1)+sum(pars) )
  # n*n/2...? or n/2...?
  return(rate)
}



Phi <- function(y,nt,t,theta,prior="swf"){
  weights_fun <- match.fun(paste(prior,"weights",sep="."))
  # prediction weights
  weights <- weights_fun(y,nt,t,theta)
  # update
  newpars <- c(NA,NA)
  newpars[1] <- theta[1]+y
  newpars[2] <- theta[2]+nt-y
  # initialise pi_dens object and iterables
  pi_dens <- rep(0,99)
  nr <- dim(weights)[1] - 1
  x1 <- seq(0.01,1-0.01,by=0.01)
  #nr <- 0:y[2]
  for (i in 0:nr){
    for (j in 0:(nt-nr)){
      if(i==0&j==0){
        pi_dens <- weights[i+1,1]*dbeta(x1,newpars[1]+i,newpars[2]+j)
      }else{
        pi_dens <- pi_dens + weights[i+1,1]*dbeta(x1,newpars[1]+i,newpars[2]+j)
      }
    }
  }
  return(list(pi_dens,newpars))
}


# input is a single distribution
# output is a mixture

# we need to extend the code to allow for a
# mixture input with corresponding mixture ouput

