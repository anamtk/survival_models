model{
  #-------------------------------------## 
  #Logistic exposure model with interval-specific covariates
  ##and custom total survival probability based on interval
  ##survival probabilities ###
  #-------------------------------------##
  
  #Ana Miller-ter Kuile and Kiona Ogle
  #March 28, 2023
  
  #Motivation for Survival Model
  # This is a survival model that has total-census
  # survival probabilities that depend on an individual living
  # or dying throughout a series of visit intervals in which status
  # was determined. Thus, while the model accounts for variation
  # in drivers of survival throughout the census period, the 
  # model relies on the distribution of the final fate data
  # for an individual (y[i]) and total census survival probability
  
  # The goal of this model is to avoid some of the ill-fitting
  # behaviors of a model that models interval level data 
  # (y[i,j]), such as is employed via the Mayfield Method 
  # (Mayfield 1975). 
  # This is likely more important when a relatively large amount
  # of the observation data are 1's (e.g. many individuals 
  # survive or individuals are surveyed many times while they're 
  # alive). This is a common data distribution for survival
  # models in ecology
  
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
  ## We calculate the unnormalized probability 
  ## of an individual surviving (q1[i]) or dying (q0[i]). Then, we 
  ## normalize this probability to get the normalized probability of 
  ## surviving by dividing q1[i]/(q1[i] + q0[i])
  #- This model requires that the dataset being imported into JAGS be sorted
  ## so that all single-interval individuals be the first set of individuals
  
  #Kelp survival specifics of the model:
  #-The model has random effects for transect nested within 
  ## site nested within the overall intercept(hierarchically centered) 
  #-The model includes a list of covariates that are dependent
  ## on kelp location, and do not change throughout the season
  #-The model includes several covariates that are dependent on
  ## the survey interval, including sst and wave power
  #-Imputed missing covariate values when missing data 
  ## are minimal
  
  #-------------------------------------## 
  # Likelihood ###
  #-------------------------------------##
  
  #nests with only one interval:
  for(i in 1:n.plants1){ #number of plants with only one interval
    
    #these data are bernoulli distributed
    #around total survey period survival
    y[i] ~ dbern(p1[i])
    
    #This overall survey period survival probability 
    #is dependent on surviving a series of previous
    #survey intervals. In this loop, we're only
    # looping through the plants that had one survey
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
  
  #for kelp with > 1 interval
  for(i in (n.plants1+1):n.plants){ 
    
    #these data are bernoulli distributed
    #around total survey period survival
    y[i] ~ dbern(p2[i])
    
    #This overall survey period survival probability 
    #is dependent on surviving a series of previous
    #survey intervals. In this loop, we're
    # looping through the nests that had more than
    # one survey interval
    
    #Unnormalized probabilities of survival and failure
    #If the plant did survive a set of intervals, it's
    # survival probability is the product of all
    # those intervals
    #y = 1
    q1[i] <- prod(p.int[i, 1:n.t[i]])
    
    #IF the plant did not survive at the end and was 
    # surveyed for more than one interval, it's survival
    # probability is the product of all the periods it did survive
    # minus the "mortality probability (1-survival) for
    # the last interval
    #y = 0, n.interval > 1
    q0[i] <- prod(p.int[i, 1:(n.t[i]-1)]) *
      (1 - p.int[i, n.t[i]])
    
    #For individuals with >1 interval, the normalized
    #survival probability is the probability of surviving
    # divided by the probability of surviving and dying
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
    
  } #multi-survey plants
  
  #Now - to get at the regression for this model,
  #we'll zoom into the daily survival probability
  # throughout the survey season and the covariates
  # that shape survival at a daily scale.
  #for all plants in dataset
  for(i in 1:n.plants){ #all plants in dataset
    for(j in 1:n.t[i]){ #interval # for plant i
      
      #each interval survival, p.int, is daily survival 
      #for that interval raised to the power of the number
      # of days in that interval
      p.int[i,j] <-  pow(ps[i,j], t[i,j])
      
      #this interval survival then goes back into 
      # the overall survival with custom
      #probabilities above in the plant loops
      
      #daily survival regression
      #daily survival = ps[i,j]
      #covariates that could be at the 
      #plant or interval scale
      logit(ps[i,j]) <- 
        #hierarchically centered random intercept for
        # transect within site
        b0.transect[Transect.num[i]] + 
        #interval covariates
        b[1]*SST[i,j] +
        b[2]*WavePower[i,j] +
        b[3]*Stipes[i,j] +
        #plant covariates
        b[4]*Diam[i] +
        b[5]*Depth[i] +
        b6[SubstrateID[i]]
    
      
    } #interval j
    
    #-------------------------------------## 
    # Imputing missing data ###
    #-------------------------------------##
    #Some covariate data are msising, so use the following to model those 
    # missing data
    #Basing these distributions off of the distributions of the 
    # data for each variable
    Depth[i] ~ dnorm(mu.d, tau.d)
    MaxDiam[i] ~dnorm(mu.md, tau.md)
    
  } #all plant
  
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  ## HIERARCHICAL EFFECTS PRIORS
  # Random covariates:
  # Hierarchical spatial random effects
  #each level depends on the level higher than it
  #Nested spatial random structure with hierarchical centering: 
  #tranescts within sites
  for(t in 1:n.transects){
    b0.transect[t] ~ dnorm(b0.site[Site.num[t]], tau.site)
  }
  
  #forests within overall intercept
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
    b6Substrate[s] ~ dnorm(0, 1E-2)
  }
  
  b6Substrate[1] <- 0
  
  #IMPUTING DATA PRIORS
  #Priors for missing covariate mean and precision
  mu.d ~ dunif(-10, 10)
  sig.d ~ dunif(0, 20)
  tau.d <- pow(sig.d, -2)
  mu.md ~ dunif(-10, 10)
  sig.md ~ dunif(0, 20)
  tau.md <- pow(sig.md, -2)

}