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
      #categorical covariates
      # Treatment categorical covariates
      b1TreatmentID[TreatmentID[i]] +
      #nest categorical covariate
      b2SpeciesID[SpeciesID[i]] +
      #continuous covariates
      b3StageID[StageID[i]] +
      #Nest continuouse covariates
      b[4]*NestHt[i]+
      b[5]*cosOrientation[i] +
      b[6]*InitDay[i]+
      #local continuouse covariates
      b[7]*Trees50[i] +
      b[8]*Trees2550[i] +
      b[9]*PercPonderosa[i] +
      #Climate covariates
      b[10]*Tmax[i] +
      b[11]*Tmax[i]^2 +
      b[12]*PPT[i] +
      b[13]*PPT[i]^2 +
      #landscape continuous covariates
      b[14]*ForestCV[i] +
      b[15]*Contag[i] +
      b[16]*OpenNm[i] +
      b[17]*LandHa[i] +
      b[18]*LandBu[i] 
    
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
    PercPonderosa[i] ~ dnorm(mu.pp, tau.pp)
    InitDay[i] ~ dnorm(mu.init, tau.init)
    cosOrientation[i] ~ dnorm(mu.orient, tau.orient)
    
    #temp is dependent on forest location
    Tmax[i] ~ dnorm(mu.tmax[Forest.num[i]], tau.tmax[Forest.num[i]])
    PPT[i] ~ dnorm(mu.ppt[Forest.num[i]], tau.ppt[Forest.num[i]])
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
  sig.transect ~ dunif(0, 50)
  sig.year ~ dunif(0, 10)
  
  tau.transect <- 1/pow(sig.transect,2)
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
  
  for(i in 4:18){
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
  
  
}