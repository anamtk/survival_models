model{
  
  for(d in 1:n.datasets){
  for(i in 1:n.indiv){
    for(j in 1:n.t[i,d]){
      
      #-------------------------------------## 
      # Likelihood ###
      #-------------------------------------##
      
      #observed values of y are of a 1/0 bernoulli distribution based on mu,
      #the period survival probability that is daily survival
      #rate raised to t, the number of days in the interval 
      
      y[i, j, d] ~ dbern(p.int[i,j, d])
      
      #period survival probability is determined from
      #regression below raised to the total number of days in the 
      # interval
      p.int[i,j,d] <- pow(ps[i,j,d], t[i,j,d])
      
      logit(ps[i,j,d]) <- 
        b0[d] + b1[d]*x[i,j,d]
      
      #-------------------------------------## 
      # Model Goodness-of-fit objects ###
      #-------------------------------------##
      
      #Create replicated data for gof
      yrep[i, j,d] ~ dbern(p.int[i,j,d])
      
      #Residuals
      resid[i,j,d] <- y[i,j,d] - p.int[i,j,d]
      
    }
    
    #to compare AUC across the board - just take the 
    #yrep and p from the final interval for each nest
    y.repkeep[i,d] <- yrep[i, n.t[i,d],d]
    p.intkeep[i,d] <- p.int[i, n.t[i,d],d]
    
    
  }
  
  #-------------------------------------## 
  # PRIORS ###
  #-------------------------------------##
  
  b0[d] ~ dnorm(0, 1E-2)
  
  b1[d] ~ dnorm(0, 1E-2)
  
  }
}