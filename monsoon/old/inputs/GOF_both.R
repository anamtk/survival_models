# Update nest survival model to have yrep and residuals
# Ana Miller-ter Kuile
# February 27, 2023

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


# Parameters to monitor ---------------------------------------------------

parms <- c("yrep", "resid")

# Load logit --------------------------------------------------------------

logit <- readRDS("/scratch/atm234/nest_survival/outputs/model_logit_MCMC_2_8_23.RDS")


# Update logit with yrep and resid ----------------------------------------

logit_update <- update(logit,
                parameters.to.save = parms,
                n.iter = 250)


saveRDS(logit_update, 
        "/scratch/atm234/nest_survival/outputs/model_logit_GOF_2_27_23.RDS")

Sys.time()
# Load cloglog ------------------------------------------------------------

cloglog <- readRDS("/scratch/atm234/nest_survival/outputs/model_cloglog_MCMC_2_8_23.RDS")


# Update cloglog with yrep and resid --------------------------------------

cloglog_update <- update(cloglog,
                         parameters.to.save = parms,
                         n.iter = 250)

saveRDS(cloglog_update, 
        "/scratch/atm234/nest_survival/outputs/model_cloglog_GOF_2_27_23.RDS")


Sys.time()

