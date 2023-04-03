#Convergence check
# Ana Miller-ter Kuile

#this script checks convergence for models run on monsoon

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

jags <- readRDS(file = '/scratch/atm234/nest_survival/outputs/model_cloglog_MCMC_2_8_23.RDS')

# Check convergence -------------------------------------------------------

mcmcplot(jags$samples,
         dir = '/scratch/atm234/nest_survival/outputs/mcmcplots/cloglog')

# Get RHat per parameter ------------------------------------------------

Rhat <- jags$Rhat

saveRDS(Rhat, '/scratch/atm234/nest_survival/outputs/cloglog_Rhat.RDS')
