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
                  "reshape2", "BayesPostEst")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

source(here("code",
            "functions",
            "GOF_functions.R"))

source(here("code",
            "functions",
            "plot_functions.R"))

# Load GOF model runs -----------------------------------------------------

mod_GOF <- readRDS(here('monsoon',
                      "empirical",
                      "model1",
                      "outputs",
                      "model1_GOF.RDS"))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- readRDS(here("data_outputs",
                      "03_JAGS_input_data",
                      "empirical",
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


# Balanced Accuracy -------------------------------------------------------

#sensitivity = true positive/(true positive + false negative)
#specificity = true negative/(true negative + false negative)
#balanced accuracy = (sensitivity+specificity)/2

y1 <- y %>%
  dplyr::select(-type) %>%
  rename("Observed_fate" = "Fate_class")

n.iter <- length(unique(yrep$Iteration))

acc_df <- data.frame(matrix(NA,
                            nrow = n.iter,
                            ncol = 3))

colnames(acc_df) <- c("sensitivity", "specificity", "accuracy")

for(i in 1:n.iter){
  acc_df[i,] <- accuracy_fun(yrep, y1, "Nest_ID", iteration.num = i)
}

(m1_acc <- acc_df %>%
  pivot_longer(1:3,
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("sensitivity", 
                                            "specificity", 
                                            "accuracy"))) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot() +
    ylim(0, 1))

