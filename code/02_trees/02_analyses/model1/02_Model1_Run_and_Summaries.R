# Get inits for survival Model 1
# Ana Miller-ter Kuile
# March 31, 2023

#this script runs an initial model for the full-survey interval 
#survival model and gets initials for monsoon runs

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  "R2jags", #jags wrapper
                  "coda",
                  "mcmcplots") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------
# 
# #load the formatted data for the JAGS model
data <- readRDS(here("data_outputs",
                     "02_trees",
                     "03_JAGS_input_data",
                     "mod1_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------


params <- c(
            'b0',
            #covariates
            'b',
            'b1'
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "02_trees",
              "02_analyses",
              "model1",
              "jags",
              "model1.R")

mod <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 4000,
                            DIC = TRUE)

mcmcplot(mod$samples)

gelman.diag(mod$samples,
            multivariate = F)

Rhat <- mod$Rhat

saveRDS(Rhat, 
        here('data_outputs',
             '02_trees',
             '04_posterior_summaries',
             'model1',
             'model1_Rhat.RDS'))

# get GOF -----------------------------------------------------------------

parms <- c("yrep", 'resid', 'p')

mod.update <- update(mod,
                     parameters.to.save = parms,
                     n.iter = 350,
                     codaOnly = TRUE)

saveRDS(mod.update, 
        here('data_outputs',
             '02_trees',
             '04_posterior_summaries',
             'model1',
             'GOF_model1.RDS'))

# Output summary stats and coda samples -----------------------------------

samples <- mod$samples

#pull all three chains out of the mcmc object
a <- as.matrix(samples[[1]])
b <- as.matrix(samples[[2]])
c <- as.matrix(samples[[3]])

#make this a dataframe bound by rows
samples_all <- as.data.frame(rbind(a,b,c))

#subset just a few for now because this is a LOT - can change the n 
# later if we want - I would probably aim for ~1000+ for any final
# analyses
samples_all2 <- slice_sample(samples_all, n = 4000)

saveRDS(samples_all2, 
        here('data_outputs',
             '02_trees',
             '04_posterior_summaries',
             'model1',
             'samples_model1.RDS'))

summary <- summary(mod$samples)

saveRDS(summary, 
        here('data_outputs',
             '02_trees',
             '04_posterior_summaries',
             'model1',
             'summary_model1.RDS'))


