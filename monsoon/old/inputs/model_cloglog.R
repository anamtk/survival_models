model {
  #-------------------------------------## 
  #Model of nest survival ###
  #-------------------------------------##
  
  for(i in 1:n.nests) { #for each nest
    for(j in 1:n.t[i]){ #and each interval in which the nest was surveyed
      
      #observed values of y are of a 1/0 bernoulli distribution based on mu,
      #the period survival probability that is daily survival
      #rate raised to t, the number of days in the interval 
      
      y[i, j] ~ dbern(mu[i,j]) 
      
      #period survival probability is determined from
      #regression below raised to the total number of days in the 
      # interval
      mu[i,j] <- pow(p[i,j], t[i,j])
      
      #daily survival probability is based on a 
      #set of covariates on probability
      cloglog(p[i, j]) <- #hierarchical structure of intercept with hierarchical
        # centering to fix identifiability issues
        b0.nest[Nest.num[i]] + #this encapsulates multiple spatial hierarchies
        #coded into the priors for this - see below
        #Crossed random effect of year:
        b0.year[Year.num[i]] + #this is summed to zero for identifiabilty
        # see in priors below
        #categorical covariates
        # Interval categorical covariate:
        b1StageID[StageID[i, j]] + 
        # Treatment categorical covariates
        b2TreatmentID[TreatmentID[i]] +
        #nest categorical covariate
        b3SpeciesID[SpeciesID[i]] +
        #continuous covariates
        #Nest continuouse covariates
        b[4]*NestHt[i]+
        b[5]*cosOrientation[i] +
        b[6]*InitDay[i]+
        #local continuouse covariates
        b[7]*Trees50[i] +
        b[8]*Trees2550[i] +
        b[9]*PercPonderosa[i] +
        #landscape continuous covariates
        b[10]*Tmax[i] +
        b[11]*Tmax[i]^2 +
        b[12]*ForestCV[i] +
        b[13]*Contag[i] +
        b[14]*OpenNm[i] +
        b[15]*LandHa[i] +
        b[16]*LandBu[i] 
      
      #-------------------------------------## 
      # Model Goodness-of-fit objects ###
      #-------------------------------------##
      
      #Create replicated data for gof
      yrep[i, j] ~ dbern(mu[i,j])
      
      #Residuals
      resid[i,j] <- y[i,j] - mu[i,j]
      
    }
    
    #-------------------------------------## 
    # Imputing missing data ###
    #-------------------------------------##
    
    #Some covariate data are msising, so use the following to model those 
    # missing data
    #Basing these distributions off of the distributions of the 
    # data for each variable
    Trees2550[i] ~ dnorm(mu.t25, tau.t25)
    Trees50[i] ~ dnorm(mu.t50, tau.t50)
    PercPonderosa[i] ~ dnorm(mu.pp, tau.pp)
    InitDay[i] ~ dnorm(mu.init, tau.init)
    cosOrientation[i] ~ dnorm(mu.orient, tau.orient)
    
    #temp is dependent on forest location
    Tmax[i] ~ dnorm(mu.tmax[Forest.num[i]], tau.tmax[Forest.num[i]])
    
  }
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  ## HIERARCHICAL EFFECTS PRIORS
  # Random covariates:
  # Hierarchical spatial random effects
  #each level depends on the level higher than it
  #Nested spatial random structure with hierarchical centering: 
  for(n in 1:n.nests){ #nests
    #nests in points effect
    b0.nest[n] ~ dnorm(b0.transect[Transect.num[n]], tau.nest)
  } 
  
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
  
  tau.nest ~ dgamma(0.1, 0.1)
  sig.nest <- 1/sqrt(tau.nest)
  
  #FIXED COVARIATE PRIORS
  #Categorical variables
  #this is all in relation to first treatment
  #Ensure treatment == 1 has the most observations!!
  for(s in 2:n.stages){
    b1StageID[s] ~ dnorm(0, 1E-2)
  }
  b1StageID[1] <- 0
  
  for(tt in 2:n.trt){
    b2TreatmentID[tt] ~ dnorm(0, 1E-2)
  }
  b2TreatmentID[1] <- 0
  
  for(s in 2:n.species){
    b3SpeciesID[s] ~ dnorm(0, 1E-2)
  }
  b3SpeciesID[1] <- 0
  
  for(i in 4:16){
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
  mu.pp ~ dunif(-10, 10)
  sig.pp ~ dunif(0, 20)
  tau.pp <- pow(sig.pp, -2)
  mu.init ~ dunif(-10, 10)
  sig.init ~ dunif(0, 20)
  tau.init <- pow(sig.init, -2)
  mu.orient ~ dunif(-10, 10)
  sig.orient ~ dunif(0, 20)
  tau.orient <- pow(sig.orient, -2)
  
  #these need to be indexed by forest ID
  for(f in 1:n.forests){
    mu.tmax[f] ~ dunif(-10, 10)
    sig.tmax[f] ~ dunif(0, 20)
    tau.tmax[f] <- pow(sig.tmax[f], -2)
  }
  
  #-------------------------------------## 
  # Covariate P-values ###
  #-------------------------------------##
  
  #generate a 1-0 vector for each covariate
  #such that 1 = + in that iteration, 0 = - in that iteration
  # the mean of this value will tell us whether something is mostly positive
  # (high mean posterior value), mostly negative (low mean posterior value)
  # or somewhree in the middle (often 0, so 0.5 mean posterior)
  
  #generates per level of categorical variables
  z.b1 <- step(b1StageID)
  z.b2 <- step(b2TreatmentID)
  z.b3 <- step(b3SpeciesID)
  
  #generate p-values for all continuous covariates
  for(i in 4:16){
    z[i] <- step(b[i])
  }
  
}

