#Monsoon script - Model3: Conditional probability interval model
# Ana Miller-ter Kuile
# April 5, 2023

#this script runs the custom interval probability model

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages
package.list <- c("jagsUI", "coda") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/survival_models/trees/model3/inputs/mod3_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(n.trees1 = data$n.trees1, 
                  n.trees = data$n.trees,
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

      

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = NULL,
                        model.file = "/scratch/atm234/survival_models/trees/model3/inputs/model3.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.burnin = 1000,
                        n.iter = 100000,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file = "/scratch/atm234/survival_models/trees/model3/outputs/model3_JAGS_model.RDS")

Sys.time()


