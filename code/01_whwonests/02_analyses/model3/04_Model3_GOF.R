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
                        "model3",
                        "outputs",
                        "model3_GOF.RDS"))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- readRDS(here("data_outputs",
                     "03_JAGS_input_data",
                     "01_whwonests",
                     "mod3_JAGS_input_data.RDS"))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  rename("Fate_class" = "data$y") %>%
  mutate(Nest_ID = 1:n(),
         type = "Observed") 

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yrep_1 <- as.data.frame(mod_GOF$sims.list$yrep_1) %>%
  mutate(Iteration = 1:n()) %>%
  pivot_longer(1:(last_col()-1),
               values_to = "yrep_1",
               names_to = "Nest_ID") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, length(Nest_ID)))

yrep_2 <- as.data.frame(mod_GOF$sims.list$yrep_2) %>%
  mutate(Iteration = 1:n()) %>%
  pivot_longer(1:(last_col()-1),
               values_to = "yrep_2",
               names_to = "Nest_ID") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, length(Nest_ID)))

yrep <- yrep_1 %>%
  full_join(yrep_2, by = c("Nest_ID", "Iteration")) %>%
  mutate(Fate_class = case_when(!is.na(yrep_1) ~ yrep_1,
                          !is.na(yrep_2) ~ yrep_2,
                          TRUE ~ NA_real_)) %>%
  dplyr::select(-yrep_1, -yrep_2) %>%
  mutate(type = "Simulated") %>%
  mutate(Nest_ID = as.numeric(Nest_ID))

# Graph observed versus simulated -----------------------------------------

#posterior predictive check graphical observation
(m3_pp <- ggplot() +
  #graph the simulated data
  geom_density(data = yrep, aes(x = Fate_class, group = Iteration, fill = type), 
               alpha = 0.2) +
  geom_density(data = y, aes(x = Fate_class, fill = type), alpha = 0.5) + 
   facet_wrap(~type))

#look like it varies by iteration how well it does



# Predicted accuracy ------------------------------------------------------


mu_p1 <- as.data.frame(mod_GOF$mean$yrep_1) %>%
  rename("P" = 'mod_GOF$mean$yrep_1')

mu_p2 <- as.data.frame(mod_GOF$mean$yrep_2) %>%
  rename("P" = 'mod_GOF$mean$yrep_2') %>%
  filter(!is.na(P))

times <- as.data.frame(data$n.t) %>%
  rename("intervals" = "data$n.t")

mu_p <- mu_p1 %>%
  bind_rows(mu_p2)

y_acc <- y %>%
  bind_cols(mu_p, times) %>%
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
80/(80+12) #87%
#1s:
224/(224+4) #98%

(mod3_acc_plot <- ggplot(y_acc, aes(x = Fate_class, y = P)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_boxplot() +
    labs(x = "Observed fate",
         y = "Predicted survival probability")  +
    annotate(geom = "text", 
             x = 0.75, y = 0.45,
             label = "87%") +
    annotate(geom = "text", 
             x = 2.25, y = 0.55,
             label = "98%"))


y_acc %>%
  group_by(Fate_class) %>%
  summarise(meanp = mean(P),
            sdp = sd(P),
            total = n(),
            sep = sdp/sqrt(total))
# AUC ---------------------------------------------------------------------


resp <- as.vector(data$y)

iteration.num <- length(mod_GOF$sims.list$p1[,1])

pred1 <- as.data.frame(t(mod_GOF$sims.list$p1))
pred2 <- as.data.frame(t(mod_GOF$sims.list$p2))
pred2 <- pred2[25:nrow(pred2),]

pred <- rbind(pred1, pred2)

mod3_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod3_AUC[i] <- AUC_JAGS3(pred, 
                          iteration.num = i, 
                          resp = resp)
}

as.data.frame(mod3_AUC) %>%
  summarise(mean = mean(mod3_AUC))

(mod3_AUC_plot <- as.data.frame(mod3_AUC) %>%
  ggplot() +
  geom_histogram(aes(x = mod3_AUC)) +
  geom_vline(xintercept = 0.96, linetype = 2) +
  labs(title = "Custom probability, logit link"))

