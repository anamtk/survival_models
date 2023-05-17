model{
  
  for(i in 1:n.indiv){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    y[i] ~ dbern(p[i])
    
    p[i] <- pow(ps[i], t[i])
    
    logit(ps[i]) <- b0 + b[1]*x1[i]
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    
    #Create replicated data for gof
    yrep[i] ~ dbern(p[i])
    
    #Residuals
    resid[i] <- y[i] - p[i]
    
    
  }
  
  #-------------------------------------## 
  # PRIORS ###
  #-------------------------------------##
  
  
  b0 ~ dnorm(0, 1E-2)
  
  for(c in 1:1){
  b[c] ~ dnorm(0, 1E-2)
  }
  
  #-------------------------------------## 
  # Covariate P-values ###
  #-------------------------------------##
  
  #generate a 1-0 vector for each covariate
  #such that 1 = + in that iteration, 0 = - in that iteration
  # the mean of this value will tell us whether something is mostly positive
  # (high mean posterior value), mostly negative (low mean posterior value)
  # or somewhree in the middle (often 0, so 0.5 mean posterior)
  

  #generate p-values for all continuous covariates
  for(i in 1:1){
    z[i] <- step(b[i])
  }
  
}