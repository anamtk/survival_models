#Monsoon script - Model1: total survey survival model
# Ana Miller-ter Kuile
# April 5, 2023

#this script runs the total survey survvial model

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
data <- readRDS("/scratch/atm234/survival_models/kelp/model1/inputs/mod1_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------

data_list <- list(#Data count variables
                  n.plants = data$n.plants,
                  n.transects = data$n.transects,
                  n.sites = data$n.sites,
                  n.substrates = data$n.substrates,
                  #Random effects IDs
                  Transect.num = data$Transect.num,
                  Site.num = data$Site.num,
                  #Covariates
                  SST = data$SST,
                  WavePower = data$WavePower,
                  Stipes = data$Stipes,
                  Diam = data$Diam,
                  Depth = data$Depth,
                  SubstrateID = data$SubstrateID,
                  #dataset
                  y = data$y,
                  #exposure
                  t = data$t)
# Parameters to save ------------------------------------------------------

params <- c(
            #Random covariate betas
            'b0.transect',
            'b0.site',
            'b0',
            #Variance/precision
            'sig.transect',
            'sig.site',
            #covariates
            'b',
            'b6'
          )

# INits -------------------------------------------------------------------

inits <- readRDS("/scratch/atm234/survival_models/kelp/model1/inputs/model1_inits.RDS")

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = inits,
                        #inits = NULL,
                        model.file = "/scratch/atm234/survival_models/kelp/model1/inputs/model1.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.burnin = 1000,
                        n.iter = 14000,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file = "/scratch/atm234/survival_models/kelp/model1/outputs/model1_JAGS_model.RDS")

Sys.time()


