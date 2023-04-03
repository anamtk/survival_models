# Assess model convergence and posteriors
# Ana Miller-ter Kuile
# November 9, 2021

#this script looks at the outputs of the nest survival model

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  'coda', #get posterior MCMC outputs
                  "bayesplot", #plot bayesian stuff
                  "tidybayes",#more bayesian plotting stuff
                  "mcmcplots", #posterior graphing
                  "patchwork",
                  "rjags")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#get model diagnostic plotting functions
source(here("code",
            "functions",
            "plot_functions.R"))

# Load model output -------------------------------------------------------

logit <- readRDS(here("monsoon",
                      "outputs",
                      "logit_Rhat.RDS"))

cloglog <- readRDS(here('monsoon',
                        'outputs',
                        'cloglog_Rhat.RDS'))
# Posterior distributions and trace plots ---------------------------------

# Gelman-Rubin ------------------------------------------------------------

parms <- c("b", "b0", "b0.nest", "b0.transect",
           "b0.year", "deviance", "sig.nest",
           "sig.transect", "sig.year",
           "b", "b1StageID", "b2TreatmentID",
           "b3TrtTime", "b4SpeciesID")

logit <- logit[which(names(logit) %in% parms)]
cloglog <- cloglog[which(names(cloglog) %in% parms)]


#both have converged
rhat_graph_fun2(logit)
rhat_graph_fun2(cloglog)




