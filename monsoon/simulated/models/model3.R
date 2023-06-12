model{
  
  #-------------------------------------## 
  # Likelihood ###
  #-------------------------------------##
  
  #individuals with only one interval
  for(d in 1:n.datasets){
  for(i in 1:n.indiv1[d]){
    
    y[i,d] ~ dbern(p1[i,d])
    
    
    #regardless of final fate (1-0), the probability of
    #surviving just one interval is just the probability
    #of survivng that interval:
    p1[i,d] <-  p.int[i,1,d]
    
    #this part of the model is equivalent to a
    # total exposure model
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_1[i,d] ~ dbern(p1[i,d])
    # 
    # #Residuals
    resid_1[i,d] <- y[i,d] - p1[i,d]
    # 
  }
  
  #for nests with >1 intervals:
  for(i in (n.indiv1[d]+1):n.indiv){
   
    y[i,d] ~ dbern(p2[i,d])
    
    #Unnormalized probabilities of survival and failure
    #If the individual did survive a set of intervals, it's
    # survival probability is the product of all
    # those intervals
    #y = 1
    q1[i,d] <- prod(p.int[i, 1:n.t[i,d],d])
    
    #IF the indiv. did not survive at the end and was 
    # surveyed for more than one interval, it's survival
    # probability is the product of all the periods it did survive
    # minus the "mortality probability (1-survival) for
    # the last interval
    #y = 0, n.interval > 1
    q0[i,d] <- prod(p.int[i, 1:(n.t[i,d]-1),d]) *
      (1 - p.int[i, n.t[i,d],d])
    
    #For individuals with >1 interval, the normalized
    #survival probability is the probability of surviving
    # divided by the probability of surviving and dying
    p2[i,d] <- q1[i,d]/(q1[i,d] + q0[i,d])
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_2[i,d] ~ dbern(p2[i,d])
    # 
    # #Residuals
    resid_2[i,d] <- y[i,d] - p2[i,d]
    # 
    
  }
  
  for(i in 1:n.indiv){
    for(j in 1:n.t[i,d]){
      
      #each interval survival, p.int, is daily survival 
      #for that interval raised to the power of the number
      # of days in that interval
      p.int[i,j,d] <-  pow(ps[i,j,d], t[i,j,d])
      
      #this interval survival then goes back into 
      # the overall nest survival with custom
      #probabilities above in the nest loops
      
      #daily survival regression
      #daily survival = ps[i,j]
      
      logit(ps[i,j,d]) <- b0[d] + b1[d]*x[i,j,d]
    }
  }
  
  b0[d] ~ dnorm(0, 1E-2)

  b1[d] ~ dnorm(0, 1E-2)
  }
  
}