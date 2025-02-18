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
                     "01_whwonests",
                     '03_JAGS_input_data',
                     "mod1_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------


params <- c(
            #Random covariate betas
            'b0.transect',
            'b0.forest',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.transect',
            'sig.forest',
            'sig.year',
            #covariates
            'b1TreatmentID',
            'b2SpeciesID',
            'b3StageID',
            'b',
            #missing data
            'sig.t25',
            'sig.t50',
            'sig.pp',
            'sig.init',
            'sig.orient',
            'sig.tmax',
            'sig.ppt'
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "01_whwonests",
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

gelman.diag(mod$samples, multivariate = F)
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
# iterations_90 iterations_95   max
# <dbl>         <dbl> <dbl>
#   1          8738        12607. 67012

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
#182

# 
# # Initials ----------------------------------------------------------------
# 
b0.forest <- mod$mean$b0.forest
sig.forest <- mod$mean$sig.forest
b0.transect <- mod$mean$b0.transect
sig.transect <- mod$mean$sig.transect
b0.year <- c(mod$mean$b0.year[1:9], NA)
sig.year <- mod$mean$sig.year
b0 <- mod$mean$b0
b1TreatmentID <- c(NA, mod$mean$b1TreatmentID[2:length(mod$mean$b1TreatmentID)])
b2SpeciesID <- c(NA, mod$mean$b2SpeciesID[2:length(mod$mean$b2SpeciesID)])
b3StageID <- c(NA, mod$mean$b3StageID[2:length(mod$mean$b3StageID)])
b <- mod$mean$b

# # Set initials ------------------------------------------------------------

inits <- list(list(b0.forest = b0.forest,
                   sig.forest = sig.forest,
                   b0.transect = b0.transect,
                   sig.transect = sig.transect,
                   b0.year = b0.year,
                   sig.year = sig.year,
                   b0 = b0,
                   b1TreatmentID = b1TreatmentID,
                   b2SpeciesID = b2SpeciesID,
                   b3StageID = b3StageID,
                   b = b),
              list(b0.forest = b0.forest + runif(length(b0.forest), min = 0, max = 1),
                   sig.forest = sig.forest + runif(length(sig.forest), min = 0, max = 1),
                   b0.transect = b0.transect +runif(length(b0.transect), min = 0, max = 1),
                   sig.transect = sig.transect +runif(length(sig.transect), min = 0, max = 1),
                   b0.year = b0.year + runif(length(b0.year), min = 0, max = 1),
                   sig.year = sig.year + runif(length(sig.year), min = 0, max = 1),
                   b0 = b0 + runif(length(b0), min = 0, max = 1),
                   b1TreatmentID =  b1TreatmentID +runif(length(b1TreatmentID), min = 0, max = 1),
                   b2SpeciesID = b2SpeciesID + runif(length(b2SpeciesID), min = 0, max = 1),
                   b3StageID = b3StageID + runif(length(b3StageID), min = 0, max = 1),
                   b = b + runif(length(b), min = 0, max = 1)),
              list(b0.forest = b0.forest - runif(length(b0.forest), min = 0, max = 1),
                   sig.forest = sig.forest + runif(length(sig.forest), min = 0, max = 1),
                   b0.transect = b0.transect -runif(length(b0.transect), min = 0, max = 1),
                   sig.transect = sig.transect +runif(length(sig.transect), min = 0, max = 1),
                   b0.year = b0.year - runif(length(b0.year), min = 0, max = 1),
                   sig.year = sig.year + runif(length(sig.year), min = 0, max = 1),
                   b0 = b0 - runif(length(b0), min = 0, max = 1),
                   b1TreatmentID =  b1TreatmentID -runif(length(b1TreatmentID), min = 0, max = 1),
                   b2SpeciesID = b2SpeciesID - runif(length(b2SpeciesID), min = 0, max = 1),
                   b3StageID = b3StageID - runif(length(b3StageID), min = 0, max = 1),
                   b = b - runif(length(b), min = 0, max = 1)))

# 
 saveRDS(inits, 
         file = here("monsoon",
                     "01_whwonests",
                     "model1",
                     "inputs",
                     "model1_inits.RDS"))
# 
# 



