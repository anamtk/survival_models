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
                        "02_kelp",
                        "model2",
                        "outputs",
                        "model2_GOF.RDS"))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- readRDS(here("data_outputs",
                     "02_kelp",
                     "03_JAGS_input_data",
                     "mod2_JAGS_input_data.RDS"))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  mutate(Plant_ID = 1:n()) %>%
  pivot_longer(1:10,
               names_to = "Interval",
               values_to = "Fate_class") %>%
  filter(!is.na(Fate_class)) %>%
  group_by(Plant_ID) %>%
  filter(Interval == max(Interval, na.rm = T)) %>%
  ungroup() %>%
  mutate(type = "Observed") 

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps <- mod_GOF$sims.list$y.repkeep

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep<- reshape2::melt(yreps) %>%
  rename("Iteration" = "Var1",
         "Plant_ID" = "Var2",
         "Fate_class" = "value") %>%
  mutate(type = "Simulated")

# Graph observed versus simulated -----------------------------------------

#posterior predictive check graphical observation
(m2_pp <- ggplot() +
   #graph the simulated data
   geom_density(data = yrep, aes(x = Fate_class, group = Iteration, fill = type), 
                alpha = 0.2) +
   geom_density(data = y, aes(x = Fate_class, fill = type), alpha = 0.5)) 


# Prediction accuracy -----------------------------------------------------

mu.p <- as.data.frame(mod_GOF$mean$p.intkeep) %>%
  rename("P" = 'mod_GOF$mean$p.intkeep')

y_acc <- y %>%
  bind_cols(mu.p) %>%
  mutate(Fate_class = as.factor(Fate_class)) %>%
  mutate(type = case_when((Fate_class == 0 & P >= .5) ~ "mis-0",
                          (Fate_class == 0 & P < .5) ~ "match-0",
                          (Fate_class == 1 & P >= 0.5) ~ "match-1",
                          (Fate_class == 1 & P < 0.5) ~ "mis-1"))


y_acc %>%
  group_by(type) %>%
  tally()

#1s:
208/(208+33)
#0s:
734/(734+863)

(mod2_acc_plot <- ggplot(y_acc, aes(x = Fate_class, y = P)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_boxplot() +
    labs(x = "Observed fate",
         y = "Predicted survival probability") +
    annotate(geom = "text", 
             x = 0.75, y = 0.45,
             label = "46%") +
    annotate(geom = "text", 
             x = 2.25, y = 0.55,
             label = "86%") )

# AUC ---------------------------------------------------------------------

resp <- as.vector(y$Fate_class)

iteration.num <- length(mod_GOF$sims.list$p.intkeep[,1])

AUC_JAGS4(mod_GOF, 
          iteration.num = 11, 
          resp = resp)

mod2_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod2_AUC[i] <- AUC_JAGS4(mod_GOF, 
                          iteration.num = i, 
                          resp = resp)
}

as.data.frame(mod2_AUC) %>%
  summarise(mean = mean(mod2_AUC))

(mod2_AUC_plot <- as.data.frame(mod2_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod2_AUC)) +
    geom_vline(xintercept = 0.74, linetype = 2) +
    labs(title = "Interval-level response \n (last survey data only), logit link"))


#THIS IS for all - would need to track pint though for all intervals
t <- mod_GOF$sims.list$p.int

layers <- dim(t)[[3]]

dfs <- lapply(1:layers,
              function(x){
                return(as.data.frame(t[,,x]))
              } )

dfs1 <- dfs %>%
  map(~mutate(., iteration = 1:n()))

full_df <- bind_rows(dfs1, .id = "interval") %>%
  pivot_longer(cols = 2:(last_col()-1),
               values_to = "p",
               names_to = "Nest_ID") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, length(Nest_ID))) %>%
  unite(col = "ID_interval",
        c("Nest_ID", "interval"),
        sep = "_") %>%
  filter(!is.na(p))


resp <- as.data.frame(data$y) %>%
  mutate(Nest_ID = 1:n()) %>%
  pivot_longer(cols = 1:(last_col()-1),
               values_to = "resp",
               names_to = "interval") %>%
  filter(!is.na(resp)) %>%
  unite(col = "ID_interval",
        c("Nest_ID", "interval"),
        sep = "_")


AUC_JAGS2(df = full_df,
          iteration.num = 3,
          resp = resp$resp)

iteration.num <- length(unique(full_df$iteration))

mod2_AUC2 <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod2_AUC2[i] <- AUC_JAGS2(df = full_df,
                           iteration.num = i,
                           resp = resp$resp)
}

as.data.frame(mod2_AUC2) %>%
  summarise(mean = mean(mod2_AUC2))

(mod2_AUC_plotall <- as.data.frame(mod2_AUC2) %>%
  ggplot() +
  geom_histogram(aes(x = mod2_AUC2)) +
  geom_vline(xintercept = 0.54, linetype = 2) +
  labs(title = "Interval-level response, logit link"))

