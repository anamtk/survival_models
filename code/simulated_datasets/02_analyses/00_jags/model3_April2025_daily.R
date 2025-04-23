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
    pi1[i,d] <-  p.int[i,1,d]
    
    #this part of the model is equivalent to a
    # total exposure model
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_1[i,d] ~ dbern(pi1[i,d])
    # 
    # #Residuals
    resid_1[i,d] <- y[i,d] - pi1[i,d]
    # 
  } #indiv1
  
  #for nests with >1 intervals:
  for(i in (n.indiv1[d]+1):n.indiv){
   
    y[i,d] ~ dbern(pi2[i,d])
    
    
    #Unnormalized probabilities of survival and failure
    #If the individual did survive a set of intervals, it's
    # survival probability is the product of all
    # those intervals
    #y = 1
    #q1[i,d] <- prod(p.int[i, 1:n.t[i,d],d])*y2[i,d]
    
    #IF the indiv. did not survive at the end and was 
    # surveyed for more than one interval, it's survival
    # probability is the product of all the periods it did survive
    # minus the "mortality probability (1-survival) for
    # the last interval
    #y = 0, n.interval > 1
    # q0[i,d] <- 1 - prod(p.int[i, 1:(n.t[i,d]-1),d])*
    #   (1-p.int[i,n.t[i,d],d])*(1-y2[i,d])
    
    #For individuals with >1 interval, the
    #probability of success is these two
    #conditions added together 
    #(first one will hold when y = 0; the other when y=1)
    pi2[i,d] <- (1 - prod(p.int[i, 1:(n.t[i,d]-1),d])*
      (1-p.int[i,n.t[i,d],d]))*(1-y2[i,d]) + 
      prod(p.int[i, 1:n.t[i,d],d])*y2[i,d]
    
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
  
  for(i in 1:n.indiv){
    for(j in 1:n.t[i,d]){

      #each interval survival, p.int, is product of
      #daily survival for all the days in that interval
      p.int[i,j,d] <- prod(ps[i,start.day[i,j,d]:end.day[i,j,d],d])

    } #n.t
    #FIX THE DAY INDEXING IT IS BROKEN GRRR
    for(k in 1:n.days[i,d]){
      #this interval survival then goes back into
      # the overall nest survival with custom
      #probabilities above in the nest loops

      #daily survival regression
      #daily survival = ps[i,j]

      logit(ps[i,k,d]) <- b0[d] + b1[d]*x[i,k,d]


      x[i,k,d] ~ dnorm(mu.x, tau.x)

    } #n.days
  } #n.indiv
    # 
    # 
    # for(i in 1:n.indiv){
    #   for(j in 1:n.t[i,d]){
    # 
    #     #each interval survival, p.int, is product of
    #     #daily survival for all the days in that interval
    #     p.int[i,j,d] <- prod(ps[i,j,1:n.days[i,j,d],d])
    # 
    #     #FIX THE DAY INDEXING IT IS BROKEN GRRR
    #     for(k in 1:n.days[i,j,d]){
    #       #this interval survival then goes back into
    #       # the overall nest survival with custom
    #       #probabilities above in the nest loops
    # 
    #       #daily survival regression
    #       #daily survival = ps[i,j]
    # 
    #       logit(ps[i,j,k,d]) <- b0[d] + b1[d]*x[i,j,k,d]
    # 
    # 
    #       x[i,j,k,d] ~ dnorm(mu.x, tau.x)
    # 
    #     } #n.days
    #   } #n.t
    # } #n.indiv
  
  b0[d] ~ dnorm(0, 1E-2)

  b1[d] ~ dnorm(0, 1E-2)
  

  
  } #ndatasets
  # 
  mu.x ~ dunif(-10, 10)
  sig.x ~ dunif(0, 20)
  tau.x <- pow(sig.x, -2)
  
  
}