model {
  #-------------------------------------## 
  #Logistic exposure model with interval-specific covariates###
  #-------------------------------------##
  
  #this script models kelp survival as the surivival 
  # of a plant dependent on covariates that describe
  #covariates that stay the same and that change throughout
  #the survey period
  
  #Attributes of the model:
  #-Data are kelp fates at each survey interval (1-0),
  ## which are Bernoulli distributed with survival probability, p
  #-The model has random effects for transect nested within
  ## overall b0 (hierarchically centered) 
  #-The model includes a list of covariates that are dependent
  ## on kelp location
  #-The model includes several covariates that are dependent on
  ## the survey interval
  #-Imputed missing covariate values when missing data 
  ## are minimal
  
  
  for(i in 1:n.plants) { #for each plant
    for(j in 1:n.t[i]){ #and each interval in which the plant was surveyed
      
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
      
      #daily survival probability is based on a 
      #set of covariates on probability
      logit(ps[i, j]) <- #hierarchical structure of intercept with hierarchical
        # centering to fix identifiability issues
        b0.transect[Transect.num[i]] + #this encapsulates multiple spatial hierarchies
        #coded into the priors for this - see below
        #Interval covariates
        b[1]*SST[i,j] +
        b[2]*WavePower[i,j] +
        b[3]*Stipes[i,j] +
        #plant covariates
        b[4]*Diam[i] +
        b[5]*Depth[i] +
        b6[SubstrateID[i]]

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
    
    #-------------------------------------## 
    # Imputing missing data ###
    #-------------------------------------##
    #Some covariate data are msising, so use the following to model those 
    # missing data
    #Basing these distributions off of the distributions of the 
    # data for each variable
    Depth[i] ~ dnorm(mu.d, tau.d)
    Diam[i] ~dnorm(mu.md, tau.md)
    
  }
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  ## HIERARCHICAL EFFECTS PRIORS
  # Random covariates:
  # Hierarchical spatial random effects
  #each level depends on the level higher than it
  #Nested spatial random structure with hierarchical centering: 
  #transects within forests
  # for(t in 1:n.transects){
  #   b0.transect[t] ~ dnorm(b0.site[Site.num[t]], tau.transect)
  # }
  # 
  # #forests within overall intercept
  # for(f in 1:n.sites){
  #   b0.site[f] ~ dnorm(b0, tau.site)
  # }
  
  #transects within overall mean
  for(t in 1:n.transects){
    b0.transect[t] ~ dnorm(b0, tau.transect)
  }

  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  #for low # of levels, from Gellman paper - define sigma
  # as uniform and then precision in relation to this sigma
  sig.transect ~ dunif(0, 10)
  #sig.site ~ dunif(0, 10)
  
  tau.transect <- 1/pow(sig.transect,2)
  #tau.site <- 1/pow(sig.site, 2)
  
  #COVARIATE PRIORS
  for(i in 1:5){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #for identifiability - define 
  #the substrate code with most obs as 0 and all others
  # as normal distributions in relation to that one
  #this is similar approach to frequentist where you set one
  #level as a reference to compare effects of others to
  #cell-reference approach:
  for(s in 2:n.substrates){
    b6[s] ~ dnorm(0, 1E-2)
  }
  
  b6[1] <- 0
  
  

  #IMPUTING DATA PRIORS
  #Priors for missing covariate mean and precision
  mu.d ~ dunif(-10, 10)
  sig.d ~ dunif(0, 20)
  tau.d <- pow(sig.d, -2)
  mu.md ~ dunif(-10, 10)
  sig.md ~ dunif(0, 20)
  tau.md <- pow(sig.md, -2)
  
}

