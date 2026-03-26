# Posterior predictive checks
# April 6, 2023
# Ana Miller-ter Kuile

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
                        'model2',
                        'GOF_model2.RDS'))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- readRDS(here("data_outputs",
                     "02_trees",
                     "03_JAGS_input_data",
                     "mod2_JAGS_input_data.RDS"))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe

#this is pulling only the last interval - is this what we want?
y <- as.data.frame(data$y) %>%
  mutate(Tree_ID = 1:n()) %>%
  pivot_longer(1:3,
               names_to = "Interval",
               values_to = "Fate_class") %>%
  filter(!is.na(Fate_class)) %>%
  group_by(Tree_ID) %>%
  filter(Interval == max(Interval, na.rm = T)) %>%
  ungroup()  

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps <- mod_GOF$sims.list$y.repkeep

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep<- reshape2::melt(yreps) %>%
  rename("sample" = "Var1",
         "Tree_ID" = "Var2",
         "Fate_class" = "value") 
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

mean <- as.data.frame(mod2_AUC) %>%
  summarise(mean = mean(mod2_AUC)) %>%
  as_vector()

(mod2_AUC_plot <- as.data.frame(mod2_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod2_AUC)) +
    geom_vline(xintercept = mean, linetype = 2) +
    labs(title = "Interval-level response \n (last survey data only), logit link"))

saveRDS(mod2_AUC,
        here('data_outputs',
             '02_trees',
             '05_for_plotting',
             'mod2AUC_last.RDS'))
# 
# #THIS IS for all - would need to track pint though for all intervals
t <- mod_GOF$sims.list$p.int
# 
layers <- dim(t)[[3]]
# 
dfs <- lapply(1:layers,
              function(x){
                return(as.data.frame(t[,,x]))
              } )
# 
dfs1 <- dfs %>%
  map(~mutate(., iteration = 1:n()))
# 
full_df <- bind_rows(dfs1, .id = "interval") %>%
  pivot_longer(cols = 2:(last_col()-1),
               values_to = "p",
               names_to = "Tree_ID") %>%
  mutate(Tree_ID = str_sub(Tree_ID, 2, length(Tree_ID))) %>%
  unite(col = "ID_interval",
        c("Tree_ID", "interval"),
        sep = "_") %>%
  filter(!is.na(p))


resp <- as.data.frame(data$y) %>%
  mutate(Tree_ID = 1:n()) %>%
  pivot_longer(cols = 1:(last_col()-1),
               values_to = "resp",
               names_to = "interval") %>%
  filter(!is.na(resp)) %>%
  unite(col = "ID_interval",
        c("Tree_ID", "interval"),
        sep = "_")

# 
AUC_JAGS2(df = full_df,
          iteration.num = 3,
          resp = resp$resp)
# 
iteration.num <- length(unique(full_df$iteration))

mod2_AUC2 <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod2_AUC2[i] <- AUC_JAGS2(df = full_df,
                           iteration.num = i,
                           resp = resp$resp)
}

mean2 <- as.data.frame(mod2_AUC2) %>%
  summarise(mean = mean(mod2_AUC2)) %>%
  as_vector()

(mod2_AUC_plotall <- as.data.frame(mod2_AUC2) %>%
  ggplot() +
  geom_histogram(aes(x = mod2_AUC2)) +
  geom_vline(xintercept = mean2, linetype = 2) +
  labs(title = "Interval-level response, logit link"))

saveRDS(mod2_AUC2,
        here('data_outputs',
             '02_trees',
             '05_for_plotting',
             'mod2AUC_all.RDS'))

# Confusion matrix --------------------------------------------------------

#last interval only
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
                           'confusion_mod2_last.RDS'))


# All intervals confusion -------------------------------------------------


#all intervals
#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps <- mod_GOF$sims.list$yrep

yrep2 <- as.data.frame.table(yreps) %>%
  rename('sample' = "Var1",
         'Tree_ID' = "Var2",
         'interval' = "Var3",
         'yrep' = "Freq") %>%
  mutate(Tree_ID = as.character(as.numeric(Tree_ID)),
         interval = as.character(as.numeric(interval)))

y2 <- resp %>%
  separate(ID_interval,
           into = c("Tree_ID", 'interval'),
           sep = "_") %>%
  rename(y = resp)

conf_df2 <- yrep2 %>%
  left_join(y2, by = c("Tree_ID", "interval"))

conf_summary2 <- conf_df2 %>%
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

ggplot(conf_summary2, aes(x = 1, y = n, color = Classification)) +
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

saveRDS(conf_summary2, here("data_outputs", 
                           '02_trees',
                           '05_for_plotting',
                           'confusion_mod2_all.RDS'))


