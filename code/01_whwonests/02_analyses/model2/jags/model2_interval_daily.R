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
 
      #each interval survival, p.int, is product
      #of daily survival for all the days in that 
      #interval
          
      p.int[i,j] <- prod(ps[i,start.day[i,j]:end.day[i,j]])
      
      #-------------------------------------## 
      # Model Goodness-of-fit objects ###
      #-------------------------------------##
      
      #Create replicated data for gof
      yrep[i, j] ~ dbern(p.int[i,j])
      
      #Residuals
      resid[i,j] <- y[i,j] - p.int[i,j]
      
    }
    
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
  sig.forest ~ dunif(0, 10)
  sig.year ~ dunif(0, 10)
  
  tau.forest <- 1/pow(sig.forest, 2)
  tau.year <- 1/pow(sig.year, 2)
  
  #FIXED COVARIATE PRIORS

  #all continuous covariate b's
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

