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
                     "03_trees",
                     "03_JAGS_input_data",
                     "mod2_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------


params <- c(
            'b0',
            #covariates
            'b',
            'b10'
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_trees",
              "02_analyses",
              "model2",
              "jags",
              "model2.R")

mod <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 4000,
                            DIC = TRUE)

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
#   1        53587.        70199. 85408

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
#225

# 
# # Initials ----------------------------------------------------------------
# 
b0 <- mod$mean$b0
b14 <- c(NA, mod$mean$b14[2:length(mod$mean$b14)])
b <- mod$mean$b

# # Set initials ------------------------------------------------------------

inits <- list(list(b0 = b0,
                   b14 = b14,
                   b = b),
              list(b0 = b0 + runif(length(b0), min = 0, max = 1),
                   b14 =  b14 +runif(length(b14), min = 0, max = 1),
                   b = b + runif(length(b), min = 0, max = 1)),
              list(b0 = b0 - runif(length(b0), min = 0, max = 1),
                   b14 =  b14 -runif(length(b14), min = 0, max = 1),
                   b = b - runif(length(b), min = 0, max = 1)))

# 
saveRDS(inits, 
        file = here("monsoon",
                    "03_trees",
                    "model2",
                    "inputs",
                    "model2_inits.RDS"))
# 
# 



