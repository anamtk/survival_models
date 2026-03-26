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

mod_GOF <- readRDS(here('data_outputs',
                        '02_trees',
                        '04_posterior_summaries',
                        'model1',
                        'GOF_model1.RDS'))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- readRDS(here("data_outputs",
                     "02_trees",
                      "03_JAGS_input_data",
                      "mod1_JAGS_input_data.RDS"))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  rename("Fate_class" = "data$y") %>%
  mutate(Tree_ID = 1:n(),
         type = "Observed") 

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps <- mod_GOF$sims.list$yrep

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep<- reshape2::melt(yreps) %>%
  rename("sample" = "Var1",
         "Tree_ID" = "Var2",
         "Fate_class" = "value") %>%
  mutate(type = "Simulated")

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

saveRDS(mod1_AUC,
        here('data_outputs',
     '02_trees',
     '05_for_plotting',
     'mod1AUC.RDS'))

# Confusion matrix --------------------------------------------------------

yrep2 <- yrep %>%
  dplyr::select(sample, Tree_ID,
                Fate_class) %>%
  rename('yrep' = "Fate_class")

y2 <- y %>%
  dplyr::select(Tree_ID,
                Fate_class) %>%
  rename('y' = "Fate_class")
  
conf_df <- yrep2 %>%
  left_join(y2, by = "Tree_ID")

conf_summary <- conf_df %>%
  mutate(category = case_when((yrep == 1 & y == 1) ~ "True Positive",
                              (yrep == 0 & y == 0) ~ "True Negative",
                              (yrep == 0 & y == 1) ~ "False Negative",
                              (yrep == 1 & y == 0) ~ "False Positive",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(category)) %>%
  group_by(sample, category) %>%
  tally() %>%
  ungroup() %>%
  mutate(Actually = case_when(category %in% c("True Positive",
                                              "False Negative") ~ "Actually Positive",
                              category %in% c("False Positive",
                                              "True Negative") ~ "Actually Negative",
                              TRUE ~ NA_character_),
         Predicted = case_when(category %in% c("True Positive",
                                               "False Positive") ~ "Predicted Positive",
                               category %in% c("False Negative",
                                               "True Negative") ~ "Predicted Negative"),
         Classification = case_when(str_detect(category, "True") ~ "Correct",
                                    str_detect(category, "False") ~ "Incorrect"))

ggplot(conf_summary, aes(x = 1, y = n, color = Classification)) +
  geom_boxplot() +
  facet_grid(Predicted ~ Actually, 
             scales = "free_x") +
  scale_y_log10() +
  scale_color_manual(values = c('#4d9221', 
                                '#c51b7d')) +
  theme(strip.background = element_blank(),
        axis.text.x  =element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  labs(y= "Count")

saveRDS(conf_summary, here("data_outputs", 
                           '02_trees',
                           '05_for_plotting',
                           'confusion_mod1.RDS'))



