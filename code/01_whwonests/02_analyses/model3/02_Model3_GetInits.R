# Running the nest survival model
# Ana Miller-ter Kuile
# November 4, 2021

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

data <- readRDS(here("data_outputs", 
                     '03_JAGS_input_data',
                     "01_whwonests",
                     "mod3_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------


params <- c(
            #Random covariate betas
            'b0.nest',
            'b0.transect',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.nest',
            'tau.nest',
            'sig.transect',
            'sig.year',
            'b',
            'b1StageID',
            'b2TreatmentID',
            'b3SpeciesID'
          )
                        



# JAGS model --------------------------------------------------------------

model <- here("code", 
              "01_whwonests",
              "02_analyses",
              "model3",
              "jags",
              "model3.R")

Sys.time()
mod <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 77000,
                            n.burnin = 1000,
                            DIC = TRUE)
Sys.time()
saveRDS(mod, here("monsoon",
                  "01_whwonests",
                  "model3",
                  "outputs",
                  "model3_JAGS_model.RDS"))
mcmcplot(mod$samples)

# Raftery -----------------------------------------------------------------

raf <- raftery.diag(mod$samples)

names <- rownames(raf[[1]]$resmatrix)
ch1 <- raf[[1]]$resmatrix[,2]
ch2 <- raf[[2]]$resmatrix[,2]
ch3 <- raf[[3]]$resmatrix[,2]

raf_all <- as.data.frame(cbind(names, 
                               ch1, ch2, ch3)) %>%
  mutate(ch1 = as.numeric(ch1),
         ch2 = as.numeric(ch2),
         ch3 = as.numeric(ch3)) %>%
  filter(!str_detect(names, "z")) %>%
  filter(!str_detect(names, "wA")) %>%
  filter(!str_detect(names, "wB")) %>%
  pivot_longer(ch1:ch3,
               names_to = "chain",
               values_to = 'iterations') 

ggplot(raf_all, aes(x = iterations/3)) +
  geom_histogram() 

raf_all %>%
  summarise(iterations_90 = quantile(iterations, 
                                     probs = 0.9, 
                                     na.rm = T)/3,
            iterations_95 = quantile(iterations,
                                     probs = 0.95,
                                     na.rm = T)/3,
            max = max(iterations, 
                      na.rm = T)/3)
# A tibble: 1 × 3
# iterations_90 iterations_95    max
# <dbl>         <dbl>  <dbl>
#   1         7339.         9634. 53885.

bu1 <- raf[[1]]$resmatrix[,1]
bu2 <- raf[[2]]$resmatrix[,1]
bu3 <- raf[[3]]$resmatrix[,1]

burn <- as.data.frame(cbind(names, bu1, bu2, bu3)) %>%
  mutate(bu1 = as.numeric(bu1),
         bu2 = as.numeric(bu2),
         bu3 = as.numeric(bu3)) %>%
  filter(!str_detect(names, "z")) %>%
  filter(!str_detect(names, "wA")) %>%
  filter(!str_detect(names, "wB")) %>%
  pivot_longer(bu1:bu3,
               names_to = "chain",
               values_to = 'iterations') 

burn %>%
  summarise(max(iterations, na.rm = T))
#643

# 
# # Initials ----------------------------------------------------------------
# 
b0.transect <- mod$mean$b0.transect
sig.transect <- mod$mean$sig.transect
b0.nest <- mod$mean$b0.nest
tau.nest <- mod$mean$tau.nest
b0.year <- c(mod$mean$b0.year[1:9], NA)
sig.year <- mod$mean$sig.year
b0 <- mod$mean$b0
b2TreatmentID <- c(NA, mod$mean$b2TreatmentID[2:length(mod$mean$b2TreatmentID)])
b1StageID <- c(NA, mod$mean$b1StageID[2:length(mod$mean$b1StageID)])
b3SpeciesID <- c(NA, mod$mean$b3SpeciesID[2:length(mod$mean$b3SpeciesID)])
b <- mod$mean$b

# Set initials ------------------------------------------------------------

inits <- list(list(b0.transect = b0.transect,
                   sig.transect = sig.transect,
                   b0.year = b0.year,
                   sig.year = sig.year,
                   b0.nest = b0.nest,
                   tau.nest = tau.nest,
                   b0 = b0,
                   b2TreatmentID = b2TreatmentID,
                   b1StageID = b1StageID,
                   b3SpeciesID = b3SpeciesID,
                   b = b),
              list(b0.transect = b0.transect +runif(length(b0.transect)),
                   sig.transect = sig.transect +runif(length(sig.transect)),
                   b0.year = b0.year + runif(length(b0.year)),
                   sig.year = sig.year + runif(length(sig.year)),
                   b0.nest = b0.nest + runif(length(b0.nest)),
                   tau.nest = tau.nest + runif(length(tau.nest)),
                   b0 = b0 + runif(length(b0)),
                   b2TreatmentID =  b2TreatmentID +runif(length(b2TreatmentID)),
                   b1StageID = b1StageID +runif(length(b1StageID)),
                   b3SpeciesID = b3SpeciesID + runif(length(b3SpeciesID)),
                   b = b + runif(length(b))),
              list(b0.transect = b0.transect -runif(length(b0.transect)),
                   sig.transect = sig.transect +runif(length(sig.transect)),
                   b0.year = b0.year - runif(length(b0.year)),
                   sig.year = sig.year + runif(length(sig.year)),
                   b0.nest = b0.nest - runif(length(b0.nest)),
                   tau.nest = tau.nest + runif(length(tau.nest)),
                   b0 = b0 - runif(length(b0)),
                   b2TreatmentID =  b2TreatmentID -runif(length(b2TreatmentID)),
                   b1StageID = b1StageID -runif(length(b1StageID)),
                   b3SpeciesID = b3SpeciesID - runif(length(b3SpeciesID)),
                   b = b - runif(length(b))))


saveRDS(inits,
        file = here("monsoon",
                    "01_whwonests",
                    "model3",
                    "inputs",
                    "model3_inits.RDS"))


