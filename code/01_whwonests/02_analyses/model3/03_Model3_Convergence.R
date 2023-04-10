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
            "00_functions",
            "plot_functions.R"))

# Load model output -------------------------------------------------------

Rhat <- readRDS(here("monsoon",
                     "01_whwonests",
                     "model3",
                     "outputs",
                     "model3_Rhat.RDS"))


# Gelman-Rubin ------------------------------------------------------------

#both have converged
rhat_graph_fun2(Rhat)
rhat_graph_fun(Rhat)




