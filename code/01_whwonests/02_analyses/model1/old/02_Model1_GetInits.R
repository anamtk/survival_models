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

#get model diagnostic plotting functions
source(here("code",
            "00_functions",
            "plot_functions.R"))

source(here("code",
            "00_functions",
            "GOF_functions.R"))

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
            'b0.forest',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.forest',
            'sig.year',
            #covariates
            'b'
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


# Graph Rhat --------------------------------------------------------------

Rhat <- mod$Rhat

#both have converged
rhat_graph_fun2(Rhat)
rhat_graph_fun(Rhat)


# Update for goodness of fit ----------------------------------------------

parms <- c("yrep", 'resid', 'p')

mod_GOF <- update(mod,
                     parameters.to.save = parms,
                     n.iter = 350,
                     codaOnly = TRUE)

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  rename("Fate_class" = "data$y") %>%
  mutate(Nest_ID = 1:n(),
         type = "Observed") 

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps <- mod_GOF$sims.list$yrep

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep<- reshape2::melt(yreps) %>%
  rename("Iteration" = "Var1",
         "Nest_ID" = "Var2",
         "Fate_class" = "value") %>%
  mutate(type = "Simulated")

# Predictive accuracy -----------------------------------------------------

mu_p <- as.data.frame(mod_GOF$mean$p) %>%
  rename("P" = 'mod_GOF$mean$p')

y_acc <- y %>%
  bind_cols(mu_p) %>%
  mutate(Fate_class = as.factor(Fate_class)) %>%
  mutate(type = case_when((Fate_class == 0 & P >= .5) ~ "mis-0",
                          (Fate_class == 0 & P < .5) ~ "match-0",
                          (Fate_class == 1 & P >= 0.5) ~ "match-1",
                          (Fate_class == 1 & P < 0.5) ~ "mis-1"))

y_acc %>%
  group_by(type) %>%
  tally()

#accuary
#0s:
zeros <- round(41/(41+43), digits = 2)*100
#1s:
ones <- round(225/(225+13), digits = 2)*100

(mod1_acc_plot <- ggplot(y_acc, aes(x = Fate_class, y = P)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_boxplot() +
    labs(x = "Observed fate",
         y = "Predicted survival probability") +
    annotate(geom = "text", 
             x = 0.75, y = 0.45,
             label = paste(zeros, "%", sep = "")) +
    annotate(geom = "text", 
             x = 2.25, y = 0.55,
             label = paste(ones, "%", sep = "")) )

y_acc %>%
  group_by(Fate_class) %>%
  summarise(meanp = mean(P),
            sdp = sd(P),
            total = n(),
            sep = sdp/sqrt(total))

# AUC ---------------------------------------------------------------------

resp <- as.vector(y$Fate_class)

iteration.num <- length(mod_GOF$sims.list$p[,1])

mod1_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod1_AUC[i] <- AUC_JAGS(mod_GOF, 
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
    labs(title = "Total survey exposure") )


# Save summaries ----------------------------------------------------------

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
             '01_whwonests',
             '04_posterior_summaries',
             'Model1_posterior_samples.RDS'))
        
summary <- summary(mod$samples)

saveRDS(summary, 
        here('data_outputs',
             '01_whwonests',
             '04_posterior_summaries',
             'Model1_posterior_summary.RDS'))



