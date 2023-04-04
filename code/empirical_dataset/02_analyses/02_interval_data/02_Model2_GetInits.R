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
                     "empirical",
                     "mod2_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #Random covariate betas
            'b0.nest',
            'b0.transect',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.nest',
            'sig.transect',
            'sig.year',
            'b',
            'b1StageID',
            'b2TreatmentID',
            'b3SpeciesID'
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "empirical_dataset",
              "02_analyses",
              "02_interval_data",
              "jags",
              "model2.R")

Sys.time()
jags <- jagsUI::jags(data = all_data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 4000,
                            DIC = TRUE)
Sys.time()
mcmcplot(LogExp.WHWO$samples)

# Raftery -----------------------------------------------------------------

raf <- raftery.diag(LogExp.WHWO$samples)

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
# iterations_90 iterations_95   max
# <dbl>         <dbl> <dbl>
#   1        14720.         22528 71705

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
#357

# 
# # Initials ----------------------------------------------------------------
# 
# b0.transect <- model_s$mean$b0.transect
# sig.transect <- model_s$mean$sig.transect
# b0.nest <- model_s$mean$b0.nest
# tau.nest <- model_s$mean$tau.nest
# b0.year <- c(model_s$mean$b0.year[1:9], NA)
# sig.year <- model_s$mean$sig.year
# b0 <- model_s$mean$b0
# b2TreatmentID <- c(NA, model_s$mean$b2TreatmentID[2:length(model_s$mean$b2TreatmentID)])
# b3TrtTime <- c(NA, model_s$mean$b3TrtTime[2:length(model_s$mean$b3TrtTime)])
# b4SpeciesID <- c(NA, model_s$mean$b4SpeciesID[2:length(model_s$mean$b4SpeciesID)])
# b <- model_s$mean$b
# 
# # Set initials ------------------------------------------------------------
# 
# inits <- list(list(b0.transect = b0.transect,
#                    sig.transect = sig.transect,
#                    b0.year = b0.year,
#                    sig.year = sig.year,
#                    b0.nest = b0.nest,
#                    tau.nest = tau.nest,
#                    b0 = b0,
#                    b2TreatmentID = b2TreatmentID, 
#                    b3TrtTime = b3TrtTime,
#                    b4SpeciesID = b4SpeciesID,
#                    b = b),
#               list(b0.transect = b0.transect +runif(length(b0.transect)),
#                    sig.transect = sig.transect +runif(length(sig.transect)),
#                    b0.year = b0.year + runif(length(b0.year)),
#                    sig.year = sig.year + runif(length(sig.year)),
#                    b0.nest = b0.nest + runif(length(b0.nest)),
#                    tau.nest = tau.nest + runif(length(tau.nest)),
#                    b0 = b0 + runif(length(b0)),
#                    b2TreatmentID =  b2TreatmentID +runif(length(b2TreatmentID)), 
#                    b3TrtTime = b3TrtTime +runif(length(b3TrtTime)),
#                    b4SpeciesID = b4SpeciesID + runif(length(b4SpeciesID)),
#                    b = b + runif(length(b))),
#               list(b0.transect = b0.transect -runif(length(b0.transect)),
#                    sig.transect = sig.transect +runif(length(sig.transect)),
#                    b0.year = b0.year - runif(length(b0.year)),
#                    sig.year = sig.year + runif(length(sig.year)),
#                    b0.nest = b0.nest - runif(length(b0.nest)),
#                    tau.nest = tau.nest + runif(length(tau.nest)),
#                    b0 = b0 - runif(length(b0)),
#                    b2TreatmentID =  b2TreatmentID -runif(length(b2TreatmentID)), 
#                    b3TrtTime = b3TrtTime -runif(length(b3TrtTime)),
#                    b4SpeciesID = b4SpeciesID - runif(length(b4SpeciesID)),
#                    b = b - runif(length(b))))
# 
# 
# saveRDS(inits, 
#         file = here("data_outputs", 
#                     "02_monsoon",
#                     "model_input_data",
#                     "nest_survival_inits.RDS"))
# 
# 

# JAGS model --------------------------------------------------------------

model <- here("code", 
              "analyses",
              "jags",
              "model_survival_feb23cloglog.R")

Sys.time()
LogExp.WHWO2 <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            #n.adapt = 10, #set to at least 1000
                            #n.burnin = 5,
                            n.iter = 4000,
                            #n.thin = 10, #higher = more thinning - takes every n.thin= iteration,
                            DIC = TRUE)

Sys.time()

# Raftery -----------------------------------------------------------------

raf <- raftery.diag(LogExp.WHWO2$samples)

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
# # A tibble: 1 × 3
# iterations_90 iterations_95   max
# <dbl>         <dbl> <dbl>
#   1         8246.        10841. 60796

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
