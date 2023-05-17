model{
  
  #-------------------------------------## 
  # Likelihood ###
  #-------------------------------------##
  
  #individuals with only one interval
  for(i in 1:n.indiv1){ 
    y[i] ~ dbern(p1[i])
    
    
    #regardless of final fate (1-0), the probability of
    #surviving just one interval is just the probability
    #of survivng that interval:
    p1[i] <-  p.int[i,1]
    
    #this part of the model is equivalent to a
    # total exposure model
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_1[i] ~ dbern(p1[i])
    # 
    # #Residuals
    resid_1[i] <- y[i] - p1[i]
    # 
  }
  
  #for nests with >1 intervals:
  for(i in (n.indiv1+1):n.indiv){
   
    y[i] ~ dbern(p2[i])
    
    #Unnormalized probabilities of survival and failure
    #If the individual did survive a set of intervals, it's
    # survival probability is the product of all
    # those intervals
    #y = 1
    q1[i] <- prod(p.int[i, 1:n.t[i]])
    
    #IF the indiv. did not survive at the end and was 
    # surveyed for more than one interval, it's survival
    # probability is the product of all the periods it did survive
    # minus the "mortality probability (1-survival) for
    # the last interval
    #y = 0, n.interval > 1
    q0[i] <- prod(p.int[i, 1:(n.t[i]-1)]) *
      (1 - p.int[i, n.t[i]])
    
    #For individuals with >1 interval, the normalized
    #survival probability is the probability of surviving
    # divided by the probability of surviving and dying
    p2[i] <- q1[i]/(q1[i] + q0[i])
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_2[i] ~ dbern(p2[i])
    # 
    # #Residuals
    resid_2[i] <- y[i] - p2[i]
    # 
    
  }
  
  for(i in 1:n.indiv){
    for(j in 1:n.t[i]){
      
      #each interval survival, p.int, is daily survival 
      #for that interval raised to the power of the number
      # of days in that interval
      p.int[i,j] <-  pow(ps[i,j], t[i,j])
      
      #this interval survival then goes back into 
      # the overall nest survival with custom
      #probabilities above in the nest loops
      
      #daily survival regression
      #daily survival = ps[i,j]
      
      logit(ps[i,j]) <- b0 + b[1]*x1[i,j]
    }
  }
  
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