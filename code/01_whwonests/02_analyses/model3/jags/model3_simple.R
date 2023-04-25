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
  for(i in (n.nests1+1):tot.nests){ 
    
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
    p_2[i] <- prod(p.int[i, 1:n.t[i]])
    
    #IF the nest did not survive at the end and was 
    # surveyed for more than one interal, it's survival
    # probability is 1 - all the periods it did survive
    # minus the "mortality probability (1-survival) for
    # the last interval
    #y = 0, n.interval > 1
    p_3[i] <-  1 - prod(p.int[i, 1:(n.t[i]-1)]) *
      (1 - p.int[i, n.t[i]])
    
    #total nest survival is based on 
    #the "if-else" of the above
    #two conditions
    #if a condition isn't met, that part
    # will become zero added to the other probabilities
    
    #only one interval and died, p1
    p2[i] <- 
      #lived through all intervals, p_2
      (y2[i]==1)*p_2[i] +
      #died in an interval after the first one , p3
      (y2[i] == 0)*p_3[i]
    
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
  for(i in 1:tot.nests){ #all nests in dataset
    for(j in 1:n.t[i]){ #interval # for nest i
      
      #daily survival = ps[i,j]
      #covariates that could be at the 
      #nest scale or at the nest + interval level
      logit(ps[i,j]) <- 
        #hierarchically centered random intercept for
        # nest within transect
        b0.transect[Transect.num[i]] + 
        #summed to zero random intercept for survey year
        b0.year[Year.num[i]] +
        b[1]*NestHt[i] 
      
      #each interval survival, p.int, is that daily
      #survival raised to the power of the number
      # of days in that interval
      p.int[i,j] <-  ps[i,j]^t[i,j]
      
      #this interval survival then goes back into 
      # the overall nest survival with custom
      #probabilities above in the nest loops
      
    } #interval j
    
    
  } #all nest loop
  
  
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
  for(i in 1:1){
    b[i] ~ dnorm(0, 1E-2)
  }
  

}