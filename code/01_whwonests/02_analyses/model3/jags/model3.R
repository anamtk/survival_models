model{
  #-------------------------------------## 
  #Logistic exposure model with daily covariates
  ##and custom total survival probability based on interval
  ##survival probabilities ###
  #-------------------------------------##
  
  #Motivation for Survival Model
  # This is a survival model that has total-census
  # survival probabilities that depend on an individual living
  # or dying throughout a series of visit intervals in which status
  # was determined. Thus, while the model accounts for variation
  # in drivers of survival throughout the census period, the 
  # model relies on the distribution of the final fate data
  # for an individual (y[i]) and total census survival probability

  
  #General attributes of the model:
  #-Data are individual fates at the end of the survey period (1-0),
  ## which are Bernoulli distributed with survival probability, p
  ## (broken into p1 and p2 in the model to account for individuals
  ## that are surveyed 1 time and >1 time)
  #- For single-survey individuals, survival probability is equivalent
  ## to a total exposure model (ps[i]^t[i])
  #- For multi-survey individuals, survival probability is dependent
  ## on fates in a series of sub-survey periods by multiplying the 
  ## survival probabilities from each of these survey intervals. 
 
  #WHWO Nest Survival specifics of the model:
  #-The model has random effects for transect nested within 
  ## forest nested within the overall intercept(hierarchically centered) 
  ## and year (sum-to-zero) that are identifiable
  #-The model includes a list of covariates that are dependent
  ## on nest location, and do not change throughout the season,
  # including habitat variables at multiple scales
  #-The model includes several covariates that are dependent on
  ## the nest and interval/day, including temperature and 
  ## development stage of the nest contents
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
    #of survivng that interval:
    p1[i] <-  p.int[i,1]
    
    #this part of the model is equivalent to a
    # total exposure model
    
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
    
    
    #For individuals with >1 interval, the propbabilyt
    #of survival is based on the probability 
    #of surviving plus the probablity of dying
    #multiplied by an indicator that is replicated
    #data of the y data (to avoid "circular"
    #errors from JAGS)
    #such that the first one will hold when y = 0; 
    #the other when y=1
    p2[i] <- (1 - prod(p.int[i,1:(n.t[i]-1)])*
                (1-p.int[i,n.t[i]]))*(1-y2[i]) + 
      prod(p.int[i,1:n.t[i]])*y2[i]
    
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
  for(i in 1:n.nests){
    for(j in 1:n.t[i]){ #total number of possible intervals
      
      #each interval survival, p.int, is product
      #of daily survival for all the days in that 
      #interval
      p.int[i,j] <- prod(ps[i,start.day[i,j]:end.day[i,j]])
      
    } #interval survivals
  
  for(k in 1:n.days[i]){ #daily data for some covariates
    
      #this interval survival then goes back into
      # the overall nest survival with custom
      #probabilities above in the nest loops
      
      #daily survival regression
      #daily survival = ps[i,k]
      #covariates that could be at the 
      #nest scale or at the nest + daily/interval level
    
      logit(ps[i,k]) <- 
        #hierarchically centered random intercept 
        b0.forest[Forest.num[i]] + 
        #summed to zero random intercept for survey year
        b0.year[Year.num[i]] +
        #Nest continuous covariates
        b[1]*NestHt[i] +
        b[2]*cosOrientation[i] +
        b[3]*InitDay[i] +
        #local-level covariates
        b[4]*Trees50[i] +
        b[5]*Trees2550[i] +
        #climate covariates
        b[6]*Tmax[i,k] +
        b[7]*Tmax[i,k]^2 +
        b[8]*PPT[i,k] +
        b[9]*PPT[i,k]^2     
      
    } #n.days
  
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
    
  } #all nest loop
  
  
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
  sig.forest ~ dunif(0, 10)
  sig.year ~ dunif(0, 10)
  
  tau.forest <- 1/pow(sig.forest,2)
  tau.year <- 1/pow(sig.year, 2)
  
  #FIXED COVARIATE PRIORS
  #for all continuous covariates b's
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