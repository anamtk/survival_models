model{
  #-------------------------------------## 
  #Model of tree survival with overall survey interval ###
  #-------------------------------------##
  
  #this script models ponderosa pine survival as the surivival 
  # of a tree dependent on covariates that describe
  #the entire survey period. It does not account for
  #covariates that might change throughout the survey
  #time
  
  #Attributes of the model:
  #-Data are final tree fates (1-0), which are Bernoulli
  ## distributed with survival probability, p
  #-The model has random effects for _____
  #-The model includes a list of covariates that are dependent
  ## on tree ID, including habitat variables and 
  ## climate variables, at multiple scales
  #-Imputed missing covariate values when missing data 
  ## are minimal

  for(i in 1:n.trees){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    #data are final tree fates (1-0) with survival
    # probability "p"
    y[i] ~ dbern(p[i])
    
    #to convert from yearly survival. 
    #accounts for different 'exposure times'
    #in this case, exposure time is measured in years
    p[i] <- pow(ps[i], t[i])
    
    #yearly survival regression
    logit(ps[i]) <-  #set of covariates on probability
      b0 +
      b[1]*AvgDBH[i] +
      b[2]*AvgBA[i] +
      b[3]*CanopyCover[i] +
      b[4]*VPD_ds[i] +
      b[5]*VPD_fa[i] +
      b[6]*VPD_ms[i] +
      b[7]*VPD_sp[i] +
      b[8]*VPD_wt[i] +
      b[9]*SWA_ds[i] +
      b[10]*SWA_fa[i] +
      b[11]*SWA_ms[i] +
      b[12]*SWA_sp[i] +
      b[13]*SWA_wt[i] +
      b14[TreatmentID[i]]
    
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

  }
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##

  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  
  #COVARIATE PRIORS
  #all other continuous covariate b's
  for(i in 1:13){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #for identifiability - define 
  #the treatment code with most obs as 0 and all others
  # as normal distributions in relation to that one
  #this is similar approach to frequentist where you set one
  #level as a reference to compare effects of others to
  #cell-reference approach:
  for(t in 2:n.trt){
    b14[t] ~ dnorm(0, 1E-2)
  }
  
  b14[1] <- 0
  
  #-------------------------------------## 
  # Covariate P-values ###
  #-------------------------------------##
  
  #generate a 1-0 vector for each covariate
  #such that 1 = + in that iteration, 0 = - in that iteration
  # the mean of this value will tell us whether something is mostly positive
  # (high mean posterior value), mostly negative (low mean posterior value)
  # or somewhree in the middle (often 0, so 0.5 mean posterior)
  
  #generates per level of categorical variables
  z.b14 <- step(b14)
 
  #generate p-values for all continuous covariates
  for(i in 1:13){
    z[i] <- step(b[i])
  }
  
  
}