model {
  #-------------------------------------## 
  #Logistic exposure model with interval-specific covariates###
  #-------------------------------------##
  
  #this script models nest survival as the surivival 
  # of a nest dependent on covariates that describe
  #covariates that stay the same and that change throughout
  #the nesting season
  
  #Attributes of the model:
  #-Data are nest fates at each survey interval (1-0),
  ## which are Bernoulli distributed with survival probability, p
  #-The model has random effects for nest nested within
  ## transect (hierarchically centered) 
  ## and year (sum-to-zero) that are identifiable
  #-The model includes a list of covariates that are dependent
  ## on nest location, including habitat variables and 
  ## climate variables, at multiple scales
  #-The model includes several covariates that are dependent on
  ## the nest and survey interval 
  #-Imputed missing covariate values when missing data 
  ## are minimal
  
  
  for(i in 1:n.nests) { #for each nest
    for(j in 1:n.t[i]){ #and each interval in which the nest was surveyed
      
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
        #Crossed random effect of year:
        b0.year[Year.num[i]] + #this is summed to zero for identifiabilty
        # see in priors below
        #Interval covariates
        # Treatment categorical covariate
        b1TreatmentID[TreatmentID[i]] +
        #nest categorical covariate
        b2SpeciesID[SpeciesID[i]] +
        #Stage categorical covariate
        b3StageID[StageID[i, j]] + 
        #continuous covariates
        #Nest continuouse covariates
        b[4]*NestHt[i] +
        b[5]*cosOrientation[i] +
        b[6]*InitDay[i] +
        #local-level covariates
        b[7]*Trees50[i] +
        b[8]*Trees2550[i] +
        b[9]*PercPonderosa[i] +
        #climate covariates
        b[10]*Tmax[i,j] +
        b[11]*Tmax[i,j]^2 +
        b[12]*PPT[i,j] +
        b[13]*PPT[i,j]^2 +
        #Landscape covariates
        b[14]*ForestCV[i] +
        b[15]*Contag[i] +
        b[16]*OpenNm[i] +
        b[17]*LandHa[i] +
        b[18]*LandBu[i] +
        b[19]*Trees50[i]*PercPonderosa[i] +
        b[20]*Trees2550[i]*PercPonderosa[i] +
        b[21]*Trees50[i]*Tmax[i,j] +
        b[22]*Trees2550[i]*Tmax[i,j]
      
      #-------------------------------------## 
      # Imputing missing data ###
      #-------------------------------------##
      
      #Some covariate data are msising, so use the following to model those 
      # missing data
      #Basing these distributions off of the distributions of the 
      # data for each variable
      
      #temp and ppt are dependent on forest location
      Tmax[i, j] ~ dnorm(mu.tmax[Forest.ID[i]], tau.tmax[Forest.ID[i]])
      PPT[i,j]~ dnorm(mu.ppt[Forest.ID[i]], tau.ppt[Forest.ID[i]])
      
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
    Trees2550[i] ~ dnorm(mu.t25, tau.t25)
    Trees50[i] ~ dnorm(mu.t50, tau.t50)
    PercPonderosa[i] ~ dnorm(mu.pp, tau.pp)
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
  #transects within forests
  for(t in 1:n.transects){
    b0.transect[t] ~ dnorm(b0.forest[Forest.num[t]], tau.transect)
  }
  
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
  sig.transect ~ dunif(0, 10)
  sig.forest ~ dunif(0, 10)
  sig.year ~ dunif(0, 10)
  
  tau.transect <- 1/pow(sig.transect,2)
  tau.forest <- 1/pow(sig.forest, 2)
  tau.year <- 1/pow(sig.year, 2)
  
  #FIXED COVARIATE PRIORS
  #Categorical variables
  #this is all in relation to first treatment
  #Ensure treatment == 1 has the most observations!!
  for(tt in 2:n.trt){
    b1TreatmentID[tt] ~ dnorm(0, 1E-2)
  }
  b1TreatmentID[1] <- 0
  
  for(s in 2:n.species){
    b2SpeciesID[s] ~ dnorm(0, 1E-2)
  }
  b2SpeciesID[1] <- 0
  
  for(s in 2:n.stages){
    b3StageID[s] ~ dnorm(0, 1E-2)
  }
  b3StageID[1] <- 0

  #all other continuous covariate b's
  for(i in 4:22){
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
     mu.ppt[f] ~ dunif(-10, 10)
     sig.ppt[f] ~ dunif(0, 20)
     tau.ppt[f] <- pow(sig.ppt[f], -2)
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
  z.b1 <- step(b1TreatmentID)
  z.b2 <- step(b2SpeciesID)
  z.b3 <- step(b3StageID)
  
  #generate p-values for all continuous covariates
  for(i in 4:22){
    z[i] <- step(b[i])
  }
  
  
}

