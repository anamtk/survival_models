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
  # (y[i,j]), such as is employed via the program MARK 
  # (Shaffer 2004, Schmidt et al. 2010). 
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
  
  #Tree survival specifics of the model:
  #-The model has random effects for ____ 
  #-The model includes a list of covariates that are dependent
  ## on tree location and the survey interval
  #-Imputed missing covariate values when missing data 
  ## are minimal
  
  #-------------------------------------## 
  # Likelihood ###
  #-------------------------------------##
  
  #trees with only one interval:
  for(i in 1:n.trees1){ #number of trees with only one interval
    
    #these data are bernoulli distributed
    #around total survey period survival
    y[i] ~ dbern(p1[i])
    
    #This overall survey period survival probability 
    #is dependent on surviving a series of previous
    #survey intervals. In this loop, we're only
    # looping through the trees that had one survey
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
    
  } #single survey trees
  
  #for trees with > 1 interval
  for(i in (n.trees1+1):n.trees){ 
    
    #these data are bernoulli distributed
    #around total survey period survival
    y[i] ~ dbern(p2[i])
    
    #This overall survey period survival probability 
    #is dependent on surviving a series of previous
    #survey intervals. In this loop, we're
    # looping through the trees that had more than
    # one survey interval
    
    #Unnormalized probabilities of survival and failure
    #If the tree did survive a set of intervals, it's
    # survival probability is the product of all
    # those intervals
    #y = 1
    q1[i] <- prod(p.int[i, 1:n.t[i]])
    
    #IF the tree did not survive at the end and was 
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
    
  } #multi-survey trees
  
  #Now - to get at the regression for this model,
  #we'll zoom into the yearly survival probability
  # throughout the survey period and the covariates
  # that shape survival at a yearly scale.
  #for all trees in dataset
  for(i in 1:n.trees){ #all trees in dataset
    for(j in 1:n.t[i]){ #interval # for tree i
      
      #to get the interval survival, p.int, we multiply
      #all the yearly survival probabilities for tree i
      #for all years in that interval, indicated by the 
      #matrices of int.start (tree x interval matrix) 
      #and int.end (tree x interval matrix)
      p.int[i,j] <- prod(p.yearly[i,int.start[i,j]:int.end[i,j]])
      
    }
  }
  
  #ANNUAL survival probabilitiy
  for(i in 1:n.trees){
    for(t in start.year[i]:end.year[i]){ #start and end year of total survey for each tree
      
      #yearly survival regression
      #covariates that could be at the 
      #tree, interval, or yearly scale
      logit(p.yearly[i,t]) <- 
        b0 +
        b1[TreatmentID[i,t]] +
        b[2]*DBH[i,t] +
        b[3]*DBH[i,t]^2 +
        b[4]*BA[i,t] +
        b[5]*CanopyCover[i,t] +
        b[6]*maxVPD[i,t] +
        b[7]*minSWA[i,t] 
      
    } #year T
  } #all tree
  
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  
  #COVARIATE PRIORS
  for(i in 2:7){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #for identifiability - define 
  #the substrate code with most obs as 0 and all others
  # as normal distributions in relation to that one
  #this is similar approach to frequentist where you set one
  #level as a reference to compare effects of others to
  #cell-reference approach:
  for(t in 2:n.trt){
    b1[t] ~ dnorm(0, 1E-2)
  }
  
  b1[1] <- 0
  
  #-------------------------------------## 
  # Covariate P-values ###
  #-------------------------------------##
  
  #generate a 1-0 vector for each covariate
  #such that 1 = + in that iteration, 0 = - in that iteration
  # the mean of this value will tell us whether something is mostly positive
  # (high mean posterior value), mostly negative (low mean posterior value)
  # or somewhree in the middle (often 0, so 0.5 mean posterior)
  
  #generates per level of categorical variables
  z.b1 <- step(b1)
  
  #generate p-values for all continuous covariates
  for(i in 2:7){
    z[i] <- step(b[i])
  }
  
  #-------------------------------------## 
  # Derived quantity of yearly survival ###
  #-------------------------------------##

  #if you wanted to know yearly survival for each year
  #regardless of tree, e.g. to know what years survival
  #was high/low, you could calculate it quickly in the model
  
  for(t in 1:n.years){ #n.years is all the years you have data
    mean.p.yearly[t] <- mean(p.yearly[,t])
  }
  
}