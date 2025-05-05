model{
  
  #-------------------------------------## 
  # Likelihood ###
  #-------------------------------------##
  
  #individuals with only one interval
  for(d in 1:n.datasets){
  for(i in 1:n.indiv1[d]){
    
    y[i,d] ~ dbern(pi1[i,d])
    
    #regardless of final fate (1-0), the probability of
    #surviving just one interval is just the probability
    #of survivng that interval:
    pi1[i,d] <-  p.int[1,d]
    
    #this part of the model is equivalent to a
    # total exposure model
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
     
    #Create replicated data for gof
    yrep_1[i,d] ~ dbern(pi1[i,d])
     
    # #Residuals
    resid_1[i,d] <- y[i,d] - pi1[i,d]
    
  } #indiv1
  
  #for nests with >1 intervals:
  for(i in (n.indiv1[d]+1):n.indiv){
   
    y[i,d] ~ dbern(pi2[i,d])
    
    #For individuals with >1 interval, the
    #probability of success is these two
    #conditions added together 
    #(first one will hold when y = 0; the other when y=1)
    pi2[i,d] <- (1 - prod(p.int[1:(n.t[i,d]-1),d])*
      (1-p.int[n.t[i,d],d]))*(1-y2[i,d]) + 
      prod(p.int[1:n.t[i,d],d])*y2[i,d]
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_2[i,d] ~ dbern(pi2[i,d])
    # 
    # #Residuals
    resid_2[i,d] <- y[i,d] - pi2[i,d]
    # 
    
  } #indiv
  
    for(j in 1:n.interval[d]){ #total number of intervals for dataset d

      #don't need i indexing here either!
      
      #each interval survival, p.int, is product of
      #daily survival for all the days in that interval
      p.int[j,d] <- prod(ps[start.day[j,d]:end.day[j,d],d])

    } #n.t
    
    for(k in 1:n.days[d]){
      #this interval survival then goes back into
      # the overall nest survival with custom
      #probabilities above in the nest loops

      #daily survival regression
      #daily survival = ps[i,j]
      
      #change indexing here - so that it 
      #is just by dataset, rather than by individual
      #since we have one nominal X for each dataset

      logit(ps[k,d]) <- b0[d] + b1[d]*x[k,d]

    } #n.days
    
  b0[d] ~ dnorm(0, 1E-2)

  b1[d] ~ dnorm(0, 1E-2)
  
  } #ndatasets

}