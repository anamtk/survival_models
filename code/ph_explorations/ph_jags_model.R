model{
  
  #xchrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://arxiv.org/pdf/2005.05952
  for(i in 1:n){ #number of indidviduals
    for(k in 1:int.obs[i]){
      cond[i,k] <- step(time[i]-a[k+1])
      #baseline hazard for individual i in time k
      HH[i,k] <- cond[i,k] * (a[k+1]-a[k]) * lambda[k] +
        (1-cond[i,k]) * (time[i]- a[k]) * lambda[k]
    }
    #cumulative harzard function
    H[i] <- sum(HH[i, 1:int.obs[i]])
  }
  
  for(i in 1:n){
    #linear predictor
    #elinpred[i] <- exp(inprod(beta[],X[i,]))
    mu[i] <- exp(b0 + b[1]*x[i])
    #log harzard function
    logHaz[i] <- log(lambda[int.obs[i]] * mu[i])
    #log survival function
    logSurv[i] <- -H[i]*mu[i]
    
    #Zeros trick for the log-likelihood
    #delta is the "observed" data - 1 if dead at end, 0 if not
    phi[i] <- 100000 - delta[i] * logHaz[i] - logSurv[i]
    zeros[i] ~ dpois(phi[i])
  }
  
  #Priors
  for(l in 1:Nbetas){
    beta[l] ~ dnorm(0, 0.001)
  }
  
  for(k in 1:m){ #number of intervals
    lambda[k] ~ dgamma(0.01, 0.01)
  }
  
}