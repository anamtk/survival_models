model{
  
  for(d in 1:n.datasets){
  for(i in 1:n.indiv){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    y[i,d] ~ dbern(p[i,d])
    
    p[i,d] <- pow(ps[i,d], t[i,d])
    
    logit(ps[i,d]) <- b0[d] + b1[d]*x[i,d]
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    
    #Create replicated data for gof
    yrep[i,d] ~ dbern(p[i,d])
    
    #Residuals
    resid[i,d] <- y[i,d] - p[i,d]
    
  }
  
  #-------------------------------------## 
  # PRIORS ###
  #-------------------------------------##
  
  b0[d] ~ dnorm(0, 1E-2)

  b1[d] ~ dnorm(0, 1E-2)

  }
}