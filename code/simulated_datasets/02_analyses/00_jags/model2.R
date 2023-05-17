model{
  
  
  for(i in 1:n.indiv){
    for(j in 1:n.t[i]){
      
      #-------------------------------------## 
      # Likelihood ###
      #-------------------------------------##
      
      #observed values of y are of a 1/0 bernoulli distribution based on mu,
      #the period survival probability that is daily survival
      #rate raised to t, the number of days in the interval 
      
      y[i, j] ~ dbern(p.int[i,j])
      
      #period survival probability is determined from
      #regression below raised to the total number of days in the 
      # interval
      p.int[i,j] <- pow(ps[i,j], t[i,j])
      
      logit(ps[i,j]) <- 
        b0 + b[1]*x1[i,j]
      
      #-------------------------------------## 
      # Model Goodness-of-fit objects ###
      #-------------------------------------##
      
      #Create replicated data for gof
      yrep[i, j] ~ dbern(p.int[i,j])
      
      #Residuals
      resid[i,j] <- y[i,j] - p.int[i,j]
      
    }
    
    #to compare AUC across the board - just take the 
    #yrep and p from the final interval for each nest
    y.repkeep[i] <- yrep[i, n.t[i]]
    p.intkeep[i] <- p.int[i, n.t[i]]
    
    
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
  
  
  for(i in 1:1){
    z[i] <- step(b[i])
  }
  
}