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
                  "mcmcplots",
                  'pROC') #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("code",
            "00_functions",
            "GOF_functions.R"))

# Load Data ---------------------------------------------------------------
# 
# #load the formatted data for the JAGS model
data <- readRDS(here("data_outputs", 
                     '03_JAGS_input_data',
                     "01_whwonests",
                     "mod1_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------


params <- c(
            #Random covariate betas
            'b0.transect',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.transect',
            'sig.year',
            'b'
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "01_whwonests",
              "02_analyses",
              "model1",
              "jags",
              "model1_simple.R")

mod <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 8000,
                            DIC = TRUE)

mcmcplot(mod$samples)

gelman.diag(mod$samples, multivariate = F)


# GOF update --------------------------------------------------------------

parms <- c("yrep", 'resid', 'p')

mod.update1 <- update(mod,
                     parameters.to.save = parms,
                     n.iter = 350,
                     codaOnly = TRUE)


y <- as.data.frame(data$y) %>%
  rename("Fate_class" = "data$y") %>%
  mutate(Nest_ID = 1:n(),
         type = "Observed") 

resp <- as.vector(y$Fate_class)

AUC_JAGS(mod_GOF = mod.update1,
         iteration.num = 11,
         resp = resp)

iteration.num <- length(mod.update1$sims.list$p[,1])

mod1_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod1_AUC[i] <- AUC_JAGS(mod.update1, 
                          iteration.num = i, 
                          resp = resp)
}

mean <- as.data.frame(mod1_AUC) %>%
  summarise(mean = mean(mod1_AUC)) %>%
  as_vector()

(mod1_AUC_plot <- as.data.frame(mod1_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod1_AUC)) +
    geom_vline(xintercept = mean, linetype = 2) +
    labs(title = "Total survey exposure, logit link"))


# Plot parameter estimates ------------------------------------------------

sum <- summary(mod$samples)
stats1 <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var= "parameter")

m1_p <- stats1 %>%
  filter(parameter %in% c("b0", "sig.transect", 
                          "sig.year", "b")) %>%
  ggplot(aes(x = parameter, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`))


b0.t1 <- stats1 %>%
  filter(str_detect(parameter, "b0.transect")) %>%
  ggplot(aes(x = parameter, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  ylim(-0.2, 11)

b0.y1 <- stats1 %>%
  filter(str_detect(parameter, "b0.year")) %>%
  ggplot(aes(x = parameter, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  ylim(-6, 2)

