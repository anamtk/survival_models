# Run low var JAGS models
# Ana Miller-ter Kuile
# May 17, 2023

#script to run models for simulated data with low predictor variation

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  "jagsUI",
                  "coda",
                  "mcmcplots")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data1 <- readRDS(here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model1low_JAGS_data.RDS"))


data2 <- readRDS(here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model2low_JAGS_data.RDS"))


data3 <- readRDS(here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model3low_JAGS_data.RDS"))


# Set up for jags models --------------------------------------------------

parms <- c("b1", "b0")

# Run model 1 -------------------------------------------------------------

model1 <- here("code",
               "simulated_datasets",
               "02_analyses",
               "00_jags",
               "model1.R")

Sys.time()
mod1 <- jagsUI::jags(data = data1,
                     inits = NULL,
                     model.file = model1,
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 4000,
                     DIC = TRUE)
Sys.time()

mcmcplot(mod1$samples)

sum1 <- summary(mod1$samples)

saveRDS(sum1, here('data_outputs',
                   'simulated',
                   '04_posterior_summaries',
                   'mod1low_summary.RDS'))

parms1 <- c( "p")

mod1update <- update(mod1,
                     parameters.to.save = parms1,
                     n.iter = 350)

saveRDS(mod1update, here('data_outputs',
                         'simulated',
                         '04_posterior_summaries',
                         'mod1low_GOFmodel.RDS'))

# Run model 2 -------------------------------------------------------------


model2 <- here("code",
               "simulated_datasets",
               "02_analyses",
               "00_jags",
               "model2.R")

Sys.time() #17 minutes
mod2 <- jagsUI::jags(data = data2,
                     inits = NULL,
                     model.file = model2,
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 1335,
                     DIC = TRUE)

Sys.time()

mcmcplot(mod2$samples)

sum2 <- summary(mod2$samples)

saveRDS(sum2, here('data_outputs',
                   'simulated',
                   '04_posterior_summaries',
                   'mod2low_summary.RDS'))

parms2 <- c("p.intkeep", "p.int")

mod2update <- update(mod2,
                     parameters.to.save = parms2,
                     n.iter = 350)

saveRDS(mod2update, here('data_outputs',
                         'simulated',
                         '04_posterior_summaries',
                         'mod2low_GOFmodel.RDS'))
# Run model 3 -------------------------------------------------------------


model3 <- here("code",
               "simulated_datasets",
               "02_analyses",
               "00_jags",
               "model3.R")

mod3 <- jagsUI::jags(data = data3,
                     inits = NULL,
                     model.file = model3,
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 4,
                     DIC = TRUE)

mcmcplot(mod3$samples)

sum3 <- summary(mod3$samples)

saveRDS(sum3, here('data_outputs',
                   'simulated',
                   '04_posterior_summaries',
                   'mod3low_summary.RDS'))

parms3 <- c("p1", "p2")

mod3update <- update(mod3,
                     parameters.to.save = parms3,
                     n.iter = 350)

saveRDS(mod3update, here('data_outputs',
                         'simulated',
                         '04_posterior_summaries',
                         'mod3low_GOFmodel.RDS'))
