model{
  #-------------------------------------## 
  #Model of kelp survival with overall survey interval ###
  #-------------------------------------##
  
  #this script models kelp survival as the surivival 
  # of a plant dependent on covariates that describe
  #the entire survey period. It does not account for
  #covariates that might change throughout the survey
  #time
  
  #Attributes of the model:
  #-Data are final kelp fates (1-0), which are Bernoulli
  ## distributed with survival probability, p
  #-The model has random effects for transect nested within
  ## site (hierarchically centered) 
  #-The model includes a list of covariates that are dependent
  ## on kelp plant, including habitat variables and 
  ## climate variables, at multiple scales
  #-Imputed missing covariate values when missing data 
  ## are minimal

  for(i in 1:n.plants){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    #data are final kelp fates (1-0) with survival
    # probability "p"
    y[i] ~ dbern(p[i])
    
    #to convert from daily survival. 
    #accounts for different 'exposure times'
    p[i] <- pow(ps[i], t[i])
    
    #daily survival regression
    logit(ps[i]) <-  #set of covariates on probability
      # centering to fix identifiability issues
      b0.transect[Transect.num[i]] + #this encapsulates transects within forest
      # see in priors below
      b[1]*SST[i] +
      b[2]*WavePower[i] +
      b[3]*Stipes[i] +
      b[4]*Diam[i] +
      b[5]*Depth[i] +
      b6[SubstrateID[i]]
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    
    #Create replicated data for gof
    yrep[i] ~ dbern(p[i])
    
    #Residuals
    resid[i] <- y[i] - p[i]
    
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
  #transects within overall site
  for(t in 1:n.transects){
    b0.transect[t] ~ dnorm(b0.site[Site.num[t]], tau.transect)
  }
  
  #sites within overall intercept
  for(f in 1:n.sites){
    b0.site[f] ~ dnorm(b0, tau.site)
  }

  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  #for low # of levels, from Gellman paper - define sigma
  # as uniform and then precision in relation to this sigma
  sig.transect ~ dunif(0, 10)
  sig.site ~ dunif(0, 10)
  
  tau.transect <- 1/pow(sig.transect,2)
  tau.site <- 1/pow(sig.site,2)
  
  #COVARIATE PRIORS
  #all other continuous covariate b's
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