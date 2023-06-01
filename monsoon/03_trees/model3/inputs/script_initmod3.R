#Monsoon script - Model3: Conditional probability interval model
# Ana Miller-ter Kuile
# April 5, 2023

#this script runs the custom interval probability model

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages
package.list <- c("jagsUI", "coda",
                  'dplyr', 'tidyr') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/survival_models/trees/model3/inputs/mod3_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(n.trees1 = data$n.trees1, 
                  n.trees = data$n.trees,
                  n.t = data$n.t,
                  n.trt = data$n.trt,
                  DBH = data$DBH,
                  BA = data$BA,
                  CanopyCover = data$CanopyCover,
                  maxVPD = data$maxVPD,
                  minSWA = data$minSWA,
                  TreatmentID = data$TreatmentID,
                  t = data$t,
                  y = data$y)

# Parameters to save ------------------------------------------------------

params <- c(
            'b0',
            #covariates
            'b',
            'b1'
          )

      

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = NULL,
                        model.file = "/scratch/atm234/survival_models/trees/model3/inputs/model3.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.iter = 4000,
                        DIC = TRUE)


# Get number of iterations for convergence --------------------------------


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
  pivot_longer(ch1:ch3,
               names_to = "chain",
               values_to = 'iterations') 

raf_all %>%
  summarise(iterations_90 = quantile(iterations, 
                                     probs = 0.9, 
                                     na.rm = T)/3,
            iterations_95 = quantile(iterations,
                                     probs = 0.95,
                                     na.rm = T)/3,
            max = max(iterations, 
                      na.rm = T)/3)


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


# 
Sys.time()


