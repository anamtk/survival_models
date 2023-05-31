model {
  #-------------------------------------## 
  #Logistic exposure model with interval-specific covariates###
  #-------------------------------------##
  
  #this script models tree survival as the surivival 
  # of a tree dependent on covariates that describe
  #covariates that stay the same and that change throughout
  #the survey period
  
  #Attributes of the model:
  #-Data are tree fates at each survey interval (1-0),
  ## which are Bernoulli distributed with survival probability, p
  #-The model has random effects for _____
  #-The model includes a list of covariates that are dependent
  ## on tree location and survey interval
  #-Imputed missing covariate values when missing data 
  ## are minimal
  
  
  for(i in 1:n.trees) { #for each tree
    for(j in 1:n.t[i]){ #and each interval in which the tree was surveyed
      
      #-------------------------------------## 
      # Likelihood ###
      #-------------------------------------##
      
      #observed values of y are of a 1/0 bernoulli distribution based on mu,
      #the period survival probability that is yearly survival
      #rate raised to t, the number of years in the interval 
      
      y[i, j] ~ dbern(p.int[i,j]) 
      
      #period survival probability is determined from
      #regression below raised to the total number of years in the 
      # interval
      p.int[i,j] <- pow(ps[i,j], t[i,j])
      
      #yearly survival probability is based on a 
      #set of covariates on probability
      logit(ps[i, j]) <- 
        b0 +
        b[1]*DBH[i,j] +
        b[2]*BA[i,j] +
        b[3]*CanopyCover[i,j] +
        b[4]*VPD_fa[i,j] +
        b[5]*VPD_ms[i,j] +
        b[6]*VPD_sp[i,j] +
        b[7]*SWA_ds[i,j] +
        b[8]*SWA_fa[i,j] +
        b[9]*SWA_wt[i,j] +
        b10[TreatmentID[i,j]]
      #-------------------------------------## 
      # Model Goodness-of-fit objects ###
      #-------------------------------------##
      
      #Create replicated data for gof
      yrep[i, j] ~ dbern(p.int[i,j])
      
      #Residuals
      resid[i,j] <- y[i,j] - p.int[i,j]
      
    }
    
    #to compare AUC across the board - just take the 
    #yrep and p from the final interval for each tree
    y.repkeep[i] <- yrep[i, n.t[i]]
    p.intkeep[i] <- p.int[i, n.t[i]]
    
    #-------------------------------------## 
    # Imputing missing data ###
    #-------------------------------------##
    #Some covariate data are msising, so use the following to model those 
    # missing data
    #Basing these distributions off of the distributions of the 
    # data for each variable
    
  }
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  
  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  
  #COVARIATE PRIORS
  for(i in 1:9){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #for identifiability - define 
  #the substrate code with most obs as 0 and all others
  # as normal distributions in relation to that one
  #this is similar approach to frequentist where you set one
  #level as a reference to compare effects of others to
  #cell-reference approach:
  for(t in 2:n.trt){
    b10[t] ~ dnorm(0, 1E-2)
  }
  
  b10[1] <- 0
  
  
  #-------------------------------------## 
  # Covariate P-values ###
  #-------------------------------------##
  
  #generate a 1-0 vector for each covariate
  #such that 1 = + in that iteration, 0 = - in that iteration
  # the mean of this value will tell us whether something is mostly positive
  # (high mean posterior value), mostly negative (low mean posterior value)
  # or somewhree in the middle (often 0, so 0.5 mean posterior)
  
  #generates per level of categorical variables
  z.b10 <- step(b10)
  
  #generate p-values for all continuous covariates
  for(i in 1:9){
    z[i] <- step(b[i])
  }
  
}

