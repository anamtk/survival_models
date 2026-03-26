model{
  #-------------------------------------## 
  #Logistic exposure model with yearly covariates
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
  #- This model requires that the dataset being imported into JAGS be sorted
  ## so that all single-interval individuals be the first set of individuals
  
  #Tree survival specifics of the model:
  #-The model includes a list of covariates that are dependent
  ## on tree location and the survey interval and year
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
    
    #For individuals with >1 interval, the
    #probability of success is these two
    #conditions added together 
    #(first one will hold when y = 0; the other when y=1)
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
    
  } #multi-survey trees
  
  #Now - to get at the regression for this model,
  #we'll zoom into the yearly survival probability
  # throughout the survey period and the covariates
  # that shape survival at a yearly scale.
  #for all trees in dataset
  for(i in 1:n.trees){ #alltrees in dataset
    for(j in 1:n.t[i]){ #interval # for tree i
      
      #UPDATE INDEXING FOR SOME OF THIS TO SIMPLIFY
      #may not need i indexing here if all
      #the trees are from the same site and 
      #have the same covariates??
      
      #each interval survival, p.int, is product of
      #yearly survival for all the years in that interval
      p.int[i,j] <-  prod(ps[i, start.year[i,j]:end.year[i,j]])
      
    }#n.t
    
    for(k in 1:n.years[i]){
      
      #this interval survival then goes back into 
      # the overall survival with custom
      #probabilities above in the tree loops
      
      #UPDATE INDEXING FOR SOME OF THIS TO SIMPLIFY
      #may not need i indexing here if all
      #the trees are from the same site and 
      #have the same covariates??
      
      #yearly survival regression
      #yearly survival = ps[i,k]
      #covariates that could be at the 
      #tree or yearly/interval scale
      logit(ps[i,k]) <- 
        b0 +
        #treatment is at the plot level
        #yearly
        b1[TreatmentID[i,k]] +
        #DBH is basd on the tree and that interval
        b[2]*DBH[i,Interval.ID[i, k]] +
        b[3]*DBH[i,Interval.ID[i, k]]^2 +
        #basal area is dependent on the plot
        #and interval
        #b[4]*BA[i,Interval.ID[i, k]] +
        #canopy cover is based on the plot
        #but we have yearly info
        b[4]*CanopyCover[i,k] +
        #vpd is at the block level and yearly
        b[5]*maxVPD[i,k] +
        #soil water is at the plot level and yearly
        b[6]*minSWA[i,k] 
      
    } #year k
    
  } #all tree
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##

  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  
  #COVARIATE PRIORS
  for(i in 2:6){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #for identifiability - define 
  #the treatment code with most obs as 0 and all others
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
  for(i in 2:6){
    z[i] <- step(b[i])
  }
  

}