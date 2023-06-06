#Monsoon script - Model2: normal daily survival model
# Ana Miller-ter Kuile
# April 5, 2023

#this script runs the original interval daily survival model

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages
package.list <- c("jagsUI", "coda",
                  'dplyr', "tidyr") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/survival_models/trees/model2/inputs/mod2_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(n.trees = data$n.trees,
                   n.t = data$n.t,
                   n.trt = data$n.trt,
                   DBH = data$DBH,
                   BA = data$BA,
                   CanopyCover = data$CanopyCover,
                   maxVPD = data$maxVPD,
                   minSWA = data$minSWA,
                   TreatmentID = data$TreatmentID,
                   t = data$t,
                   y = data$y)

# Parameters to save ------------------------------------------------------
params <- c(
            'b0',
            #covariates
            'b',
            'b1'
          )
# Inits -------------------------------------------------------------------

inits <- readRDS("/scratch/atm234/survival_models/trees/model2/outputs/mod2_JAGS_inits.RDS")

# JAGS model --------------------------------------------------------------



mod <- jagsUI::jags(data = data_list,
                        inits = inits,
                        model.file = "/scratch/atm234/survival_models/trees/model2/inputs/model2.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.iter = 171900,
                        n.burnin = 1000,
                        DIC = TRUE)


# 
#save as an R data object
saveRDS(mod, 
        file = "/scratch/atm234/survival_models/trees/model2/outputs/model2_JAGS_model.RDS")

# 
Sys.time()


