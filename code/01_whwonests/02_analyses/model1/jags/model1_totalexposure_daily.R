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
    
    #to convert from daily survival. 
    #accounts for different 'exposure times'
    p[i] <- prod(ps[i,1:(t[i])])
    #p[i] <- pow(ps[i], t[i])
    #prod(p.int[i,1:(n.t[i]-1)])
    
    for(y in 1:t[i]){
    #daily survival regression
    logit(ps[i,y]) <-  #set of covariates on probability
      # centering to fix identifiability issues
      b0.forest[Forest.num[i]] + #this encapsulates nests within forest
      #Crossed random effect of year:
      b0.year[Year.num[i]] + #this is summed to zero for identifiabilty
      # see in priors below
      #Nest continuouse covariates
      b[1]*NestHt[i]+
      b[2]*cosOrientation[i] +
      b[3]*InitDay[i]+
      #local continuouse covariates
      b[4]*Trees50[i] +
      b[5]*Trees2550[i] +
      #Climate covariates
      b[6]*Tmax[i,y] +
      b[7]*Tmax[i,y]^2 +
      b[8]*PPT[i,y] +
      b[9]*PPT[i,y]^2
    
    }
    
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
    Trees2550[i] ~ dnorm(mu.t25, tau.t25)
    Trees50[i] ~ dnorm(mu.t50, tau.t50)
    InitDay[i] ~ dnorm(mu.init, tau.init)
    cosOrientation[i] ~ dnorm(mu.orient, tau.orient)
  }
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  ## HIERARCHICAL EFFECTS PRIORS
  # Random covariates:
  # Hierarchical spatial random effects
  #each level depends on the level higher than it
  #Nested spatial random structure with hierarchical centering: 
  #forests within overall intercept
  for(f in 1:n.forests){
    b0.forest[f] ~ dnorm(b0, tau.forest)
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
  sig.forest ~ dunif(0, 50)
  sig.year ~ dunif(0, 50)
  
  tau.forest <- 1/pow(sig.forest,2)
  tau.year <- 1/pow(sig.year, 2)
  
  #FIXED COVARIATE PRIORS
  #all  continuous covariate b's
  for(i in 1:9){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #MISSING DATA IMPUTING PRIORS
  #Priors for mean and tau of missing covariates in the model
  mu.t25 ~ dunif(-10, 10)
  sig.t25 ~ dunif(0, 20)
  tau.t25 <- pow(sig.t25, -2)
  mu.t50 ~ dunif(-10, 10)
  sig.t50 ~ dunif(0, 20)
  tau.t50 <- pow(sig.t50, -2)
  mu.init ~ dunif(-10, 10)
  sig.init ~ dunif(0, 20)
  tau.init <- pow(sig.init, -2)
  mu.orient ~ dunif(-10, 10)
  sig.orient ~ dunif(0, 20)
  tau.orient <- pow(sig.orient, -2)
  
  #-------------------------------------## 
  # Covariate P-values ###
  #-------------------------------------##
  
  #generate a 1-0 vector for each covariate
  #such that 1 = + in that iteration, 0 = - in that iteration
  # the mean of this value will tell us whether something is mostly positive
  # (high mean posterior value), mostly negative (low mean posterior value)
  # or somewhree in the middle (often 0, so 0.5 mean posterior)
  
  #generate p-values for all continuous covariates
  for(i in 1:9){
    z[i] <- step(b[i])
  }
  
  
}