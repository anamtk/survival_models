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
                        'model3',
                        'GOF_model3.RDS'))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- readRDS(here("data_outputs",
                     "02_trees",
                     "03_JAGS_input_data",
                     "mod3_JAGS_input_data.RDS"))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  rename("Fate_class" = "data$y") %>%
  mutate(Tree_ID = 1:n(),
         type = "Observed") 

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yrep_1 <- as.data.frame(mod_GOF$sims.list$yrep_1) %>%
  mutate(sample = 1:n()) %>%
  pivot_longer(1:(last_col()-1),
               values_to = "yrep_1",
               names_to = "Tree_ID") %>%
  mutate(Tree_ID = str_sub(Tree_ID, 2, length(Tree_ID)))

yrep_2 <- as.data.frame(mod_GOF$sims.list$yrep_2) %>%
  mutate(sample = 1:n()) %>%
  pivot_longer(1:(last_col()-1),
               values_to = "yrep_2",
               names_to = "Tree_ID") %>%
  mutate(Tree_ID = str_sub(Tree_ID, 2, length(Tree_ID)))

yrep <- yrep_1 %>%
  full_join(yrep_2, by = c("Tree_ID", "sample")) %>%
  mutate(Fate_class = case_when(!is.na(yrep_1) ~ yrep_1,
                          !is.na(yrep_2) ~ yrep_2,
                          TRUE ~ NA_real_)) %>%
  dplyr::select(-yrep_1, -yrep_2) %>%
  mutate(type = "Simulated") %>%
  mutate(Tree_ID = as.numeric(Tree_ID))

# AUC ---------------------------------------------------------------------

resp <- as.vector(data$y)

iteration.num <- length(mod_GOF$sims.list$p1[,1])

pred1 <- as.data.frame(t(mod_GOF$sims.list$p1))
pred2 <- as.data.frame(t(mod_GOF$sims.list$p2))
pred2 <- pred2[186:nrow(pred2),]

pred <- rbind(pred1, pred2)

AUC_JAGS3(pred, 
          iteration.num = 3, 
          resp = resp)

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

saveRDS(mod3_AUC,
        here('data_outputs',
             '02_trees',
             '05_for_plotting',
             'mod3AUC.RDS'))

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
                           'confusion_mod3.RDS'))



