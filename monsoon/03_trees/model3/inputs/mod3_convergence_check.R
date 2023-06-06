
# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "dplyr",
                  "tidyr", "ggplot2", 
                  'mcmcplots', "stringr",
                  "coda", "htmltools") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Import model ------------------------------------------------------------

mod <- readRDS(file = "/scratch/atm234/survival_models/trees/model3/outputs/model3_JAGS_model.RDS")


# Check convergence -------------------------------------------------------

mcmcplot(mod$samples,
         dir = '/scratch/atm234/survival_models/trees/model3/outputs/mcmcplots/')

# Get RHat per parameter ------------------------------------------------

Rhat <- mod$Rhat

saveRDS(Rhat, '/scratch/atm234/survival_models/trees/model3/outputs/model3_Rhat.RDS')
