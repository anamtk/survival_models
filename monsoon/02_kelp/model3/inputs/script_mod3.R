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
data <- readRDS("/scratch/atm234/survival_models/kelp/model3/inputs/mod3_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(n.plants1 = data$n.plants1, 
                  n.plants = data$n.plants,
                  n.t = data$n.t,
                  n.transects =data$n.transects,
                  n.sites = data$n.sites,
                  n.substrates = data$n.substrates,
                  Transect.num = data$Transect.num,
                  Site.num = data$Site.num,
                  SST = data$SST,
                  WavePower = data$WavePower,
                  Stipes = data$Stipes,
                  Diam = data$Diam,
                  Depth = data$Depth,
                  SubstrateID = data$SubstrateID,
                  t = data$t,
                  y = data$y)

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

      

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = NULL,
                        model.file = "/scratch/atm234/survival_models/kelp/model3/inputs/model3.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.burnin = 1000,
                        n.iter = 77300,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file = "/scratch/atm234/survival_models/kelp/model3/outputs/model3_JAGS_model.RDS")

Sys.time()


