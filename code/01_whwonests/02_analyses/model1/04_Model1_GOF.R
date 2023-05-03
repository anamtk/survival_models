# Posterior predictive checks
# April 6, 2023
# Ana Miller-ter Kuile

# this is a hack of the bayesplot functionality to generate posterior
# predictive check graphs - specifically to assess - is the model family and link
# function I've selected appropriate for the data I have, or do I need to consider
# a different link or distribution (e.g. logit vs. cloglog link for binomial data; 
# poisson vs. negative binomial distribution for count data)
#then - i generate balanced accuracy values for each iteration of the model
#to see how good the model is prediction overall (balanced accuracy) and at 
# 1s (sensitivity) and 0's (specificity) separately

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "coda", "bayesplot",
                  "jagsUI",
                  "reshape2", "BayesPostEst",
                  "pROC")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

source(here("code",
            "00_functions",
            "GOF_functions.R"))

source(here("code",
            "00_functions",
            "plot_functions.R"))

# Load GOF model runs -----------------------------------------------------

mod_GOF <- readRDS(here('monsoon',
                      "01_whwonests",
                      "model1",
                      "outputs",
                      "model1_GOF.RDS"))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- readRDS(here("data_outputs",
                     "01_whwonests",
                      "03_JAGS_input_data",
                      "mod1_JAGS_input_data.RDS"))

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

# Graph observed versus simulated -----------------------------------------

#posterior predictive check graphical observation
(m1_pp <- ggplot() +
  #graph the simulated data
  geom_density(data = yrep, aes(x = Fate_class, group = Iteration, fill = type), 
               alpha = 0.2) +
  geom_density(data = y, aes(x = Fate_class, fill = type), alpha = 0.5))

#look pretty good except the wonky iterations with higher in both sizes? so
# confusing...


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
79/(79+13)
#1s:
225/(225+3)

(mod1_acc_plot <- ggplot(y_acc, aes(x = Fate_class, y = P)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
  geom_boxplot() +
  labs(x = "Observed fate",
       y = "Predicted survival probability") +
  annotate(geom = "text", 
           x = 0.75, y = 0.45,
           label = "86%") +
  annotate(geom = "text", 
           x = 2.25, y = 0.55,
           label = "99%") )
  
y_acc %>%
  group_by(Fate_class) %>%
  summarise(meanp = mean(P),
            sdp = sd(P),
            total = n(),
            sep = sdp/sqrt(total))

# AUC ---------------------------------------------------------------------

resp <- as.vector(y$Fate_class)

AUC_JAGS(mod_GOF = mod_GOF,
         iteration.num = 11,
         resp = resp)

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
  labs(title = "Total survey exposure, logit link") )

