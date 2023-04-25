model{
  #-------------------------------------## 
  #Model of nest survival with overall survey interval ###
  #-------------------------------------##
  
  #this script models nest survival as the surivival 
  # of a nest dependent on covariates that describe
  #the entire nesting season. It does not account for
  #covariates that might change throughout the nesting
  # season
  
  #Attributes of the model:
  #-Data are final nest fates (1-0), which are Bernoulli
  ## distributed with survival probability, p
  #-The model has random effects for transect (hierarchically
  ## centered) and year (sum-to-zero) that are identifiable
  #-The model includes a list of covariates that are dependent
  ## on nest location, including habitat variables and 
  ## climate variables, at multiple scales
  #-Imputed missing covariate values when missing data 
  ## are minimal

  for(i in 1:n.nests){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    #data are final nest fates (1-0) with survival
    # probability "p"
    y[i] ~ dbern(p[i])
    
    #to convert to daily survival. 
    #accounts for different 'exposure times'
    p[i] <- ps[i]^t[i]
    
    logit(ps[i]) <-  #set of covariates on probability
      # centering to fix identifiability issues
      b0.transect[Transect.num[i]] + #this encapsulates nests on a transect
      #Crossed random effect of year:
      b0.year[Year.num[i]] + #this is summed to zero for identifiabilty
      # see in priors below
      #Nest continuouse covariates
      b[1]*NestHt[i]
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    
    #Create replicated data for gof
    yrep[i] ~ dbern(p[i])
    
    #Residuals
    resid[i] <- y[i] - p[i]

  }
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  ## HIERARCHICAL EFFECTS PRIORS
  # Random covariates:
  # Hierarchical spatial random effects
  #each level depends on the level higher than it
  #Nested spatial random structure with hierarchical centering: 
  for(t in 1:n.transects){
    b0.transect[t] ~ dnorm(b0, tau.transect)
  }
  
  #Crossed effect for year
  #this effect is summed to zero for identifiability issues
  
  #for every year but the last one:
  for(y in 1:(n.years-1)){
    b0.year[y] ~ dnorm( 0, tau.year)
  }
  #set the last year to be the -sum of all other years so the 
  # overall fo all year levels == 0
  b0.year[n.years] <- -sum(b0.year[1:(n.years-1)])
  
  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  #for low # of levels, from Gellman paper - define sigma
  # as uniform and then precision in relation to this sigma
  sig.transect ~ dunif(0, 10)
  sig.year ~ dunif(0, 10)
  
  tau.transect <- 1/pow(sig.transect,2)
  tau.year <- 1/pow(sig.year, 2)
  
  #FIXED COVARIATE PRIORS
  #Categorical variables
  #this is all in relation to first treatment
  #Ensure treatment == 1 has the most observations!!
  for(i in 1:1){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  
  
}