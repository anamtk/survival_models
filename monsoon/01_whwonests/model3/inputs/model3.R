model{
  #-------------------------------------## 
  #Logistic exposure model with interval-specific covariates
  ##and custom total survival probability based on interval
  ##survival probabilities ###
  #-------------------------------------##
  
  #Ana Miller-ter Kuile and Kiona Ogle
  #March 28, 2023
  
  # This script outlines a nest survival model
  # that has total nest survivals contingent on
  # the fact that a nest lived or died throughout
  # the nesting period.
  # The goal of this model is to avoid some of the
  # ill-fitting behaviours of the model that 
  # models interval level data (y[i,j]) when
  # very many of the y's are 1's while also allowing 
  # there to be covariates that vary (eg. climate)
  # throughout different visit intervals of the 
  # nesting season
  
  #Attributes of the model:
  #-Data are nest fates at the end of the survey period (1-0),
  ## which are Bernoulli distributed with survival probability, p
  #- The model has conditional probabilities based on final nest
  ## fate and the number of survey intervals that nest was observed in
  #-These conditional probabilities are broken up between nests that
  ## were only surveyed once and those that were surveyed > 1 time 
  ## due since the conditional probability that y==0 and intervals>1
  ## will break the code since it requires n.intervals[i]-1 to be >= 1
  ## This requires that the dataset being imported into JAGS be sorted
  ## so that all single-interval nests be the first set of nests
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
  
  #-------------------------------------## 
  # Likelihood ###
  #-------------------------------------##
  
  #nests with only one interval:
  for(i in 1:n.nests1){ #number of nests with only one interval
    
    #these data are bernoulli distributed
    #around total nesting period survival
    y[i] ~ dbern(p1[i])
    
    #This overall nesting period survival probability 
    #is dependent on surviving a series of previous
    #survey intervals. In this loop, we're only
    # looping through the nests that had one survey
    # interval
    
    #regardless of final fate (1-0), the probability of
    #surviving just one interval is just the probability
    #of survivng that interval
    p1[i] <-  p.int[i,1]
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_1[i] ~ dbern(p1[i])
    # 
    # #Residuals
    resid_1[i] <- y[i] - p1[i]
    # 
    
  } #single survey nests
  
  #for nests with > 1 interval
  for(i in (n.nests1+1):n.nests){ 
    
    #these data are bernoulli distributed
    #around total nesting period survival
    y[i] ~ dbern(p2[i])
    
    #This overall nesting period survival probability 
    #is dependent on surviving a series of previous
    #survey intervals. In this loop, we're
    # looping through the nests that had more than
    # one survey interval
    
    #If the nest did survive a set of intervals, it's
    # survival probability is the product of all
    # those intervals
    #y = 1
    q1[i] <- prod(p.int[i, 1:n.t[i]])
    
    #IF the nest did not survive at the end and was 
    # surveyed for more than one interal, it's survival
    # probability is 1 - all the periods it did survive
    # minus the "mortality probability (1-survival) for
    # the last interval
    #y = 0, n.interval > 1
    # q0[i] <-  1 - prod(p.int[i, 1:(n.t[i]-1)]) *
    #   (1 - p.int[i, n.t[i]])
    
    q0[i] <- prod(p.int[i, 1:(n.t[i]-1)]) *
      (1 - p.int[i, n.t[i]])
    
    #total nest survival is based on 
    #the "if-else" of the above
    #two conditions
    #if a condition isn't met, that part
    # will become zero added to the other probabilities
    
    #only one interval and died, p1
    # p2[i] <- 
    #   #lived through all intervals, p_2
    #   (y2[i]==1)*pa[i] +
    #   #died in an interval after the first one , p3
    #   (y2[i] == 0)*pb[i]
    
    p2[i] <- q1[i]/(q1[i] + q0[i])
    
    #-------------------------------------## 
    # Model Goodness-of-fit objects ###
    #-------------------------------------##
    # 
    #Create replicated data for gof
    yrep_2[i] ~ dbern(p2[i])
    # 
    # #Residuals
    resid_2[i] <- y[i] - p2[i]
    # 
    
  } #multi-survey nests
  
  #Now - to get at the regression for this model,
  #we'll zoom into the daily survival probability
  # throughout the nesting season and the covariates
  # that shape survival at a daily scale.
  #for all nests in dataset
  for(i in 1:n.nests){ #all nests in dataset
    for(j in 1:n.t[i]){ #interval # for nest i
      
      #each interval survival, p.int, is that daily
      #survival raised to the power of the number
      # of days in that interval
      p.int[i,j] <-  pow(ps[i,j], t[i,j])
      
      #this interval survival then goes back into 
      # the overall nest survival with custom
      #probabilities above in the nest loops
      
      #daily survival regression
      #daily survival = ps[i,j]
      #covariates that could be at the 
      #nest scale or at the nest + interval level
      logit(ps[i,j]) <- 
        #hierarchically centered random intercept for
        # transect within forest
        b0.transect[Transect.num[i]] + 
        #summed to zero random intercept for survey year
        b0.year[Year.num[i]] +
        #Treatment category at nest site
        b1TreatmentID[TreatmentID[i]] +
        #nest species categorical covariate
        b2SpeciesID[SpeciesID[i]] +
        #category of the stage the nest was in in 
        # interval j
        b3StageID[StageID[i,j]] +
        #Nest continuous covariates
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
        b[18]*LandBu[i]
      
      #-------------------------------------## 
      # Imputing missing data ###
      #-------------------------------------##
      
      #Some covariate data are msising, so use the following to model those 
      # missing data
      #Basing these distributions off of the distributions of the 
      # data for each variable
      
      #temp is dependent on forest location
      Tmax[i,j] ~ dnorm(mu.tmax[Forest.ID[i]], tau.tmax[Forest.ID[i]])
      PPT[i,j] ~ dnorm(mu.ppt[Forest.ID[i]], tau.ppt[Forest.ID[i]])
      
      
    } #interval j
    
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
    
  } #all nest loop
  
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  ## HIERARCHICAL EFFECTS PRIORS
  # Random covariates:
  # Hierarchical spatial random effects
  #each level depends on the level higher than it
  #Nested spatial random structure with hierarchical centering: 
  #tranescts within forests
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
  tau.forest <- 1/pow(sig.forest,2)
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

  #for all other continuous covariates b's
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