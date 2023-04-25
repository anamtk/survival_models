# Running the nest survival model
# Ana Miller-ter Kuile
# November 4, 2021

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

data <- readRDS(here("data_outputs", 
                     '03_JAGS_input_data',
                     "01_whwonests",
                     "mod3_JAGS_input_data.RDS"))

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
              "model3",
              "jags",
              "model3_simple.R")

Sys.time()
mod3 <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 30000,
                            DIC = TRUE)
Sys.time()

#mcmcplot(mod3$samples)
gelman.diag(mod3$samples, multivariate = F)

# GOF UPdate --------------------------------------------------------------

parms <- c("yrep_1", "yrep_2", 'resid_1', "resid_2", 'p1', "p2")

mod.update3 <- update(mod3,
                     parameters.to.save = parms,
                     n.iter = 350,
                     codaOnly = TRUE)


# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  rename("Fate_class" = "data$y") %>%
  mutate(Nest_ID = 1:n(),
         type = "Observed") 

# AUC ---------------------------------------------------------------------


resp <- as.vector(data$y)

iteration.num <- length(mod.update3$sims.list$p1[,1])

#EDITING FUNCTION HERE

pred1 <- as.data.frame(t(mod.update3$sims.list$p1))
pred2 <- as.data.frame(t(mod.update3$sims.list$p2))
pred2 <- pred2[25:nrow(pred2),]

pred <- rbind(pred1, pred2)

mod3_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod3_AUC[i] <- AUC_JAGS3(pred, 
                           iteration.num = i, 
                           resp = resp)
}

mean <- as.data.frame(mod3_AUC) %>%
  summarise(mean = mean(mod3_AUC)) %>%
  as_vector()

(mod3_AUC_plot <- as.data.frame(mod3_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod3_AUC)) +
    geom_vline(xintercept = mean, linetype = 2) +
    labs(title = "Custom probability, logit link"))


# Parameter values --------------------------------------------------------

sum3 <- summary(mod3$samples)
stats3 <- as.data.frame(sum3$quantiles) %>%
  rownames_to_column(var= "parameter")

m3_p <- stats3 %>%
  filter(parameter %in% c("b0", "sig.transect", 
                          "sig.year", "b")) %>%
  ggplot(aes(x = parameter, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`))

b0.t3 <- stats3 %>%
  filter(str_detect(parameter, "b0.transect")) %>%
  ggplot(aes(x = parameter, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`))  +
  ylim(-0.2, 11)

b0.y3 <- stats3 %>%
  filter(str_detect(parameter, "b0.year")) %>%
  ggplot(aes(x = parameter, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  ylim(-6, 2)

