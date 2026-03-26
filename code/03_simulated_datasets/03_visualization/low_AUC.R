

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'patchwork',
                  'pROC')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("code",
            "00_functions",
            "GOF_functions.R"))

theme_set(theme_bw())

#import all GOF datasets
#import all OG y datasets as well

# Y datasets --------------------------------------------------------------

data1 <- readRDS(here("data_outputs",
                      "03_simulated",
                      "03_jags_input_data",
                      "Model1low_JAGS_data.RDS"))

y1low <- as.data.frame(data1$y)%>%
  rownames_to_column(var = "individual") %>%
  pivot_longer(2:last_col(),
               names_to = "site", 
               values_to = "y") %>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual))

data2 <- readRDS(here("data_outputs",
                           "03_simulated",
                           "03_jags_input_data",
                           "Model2low_JAGS_data.RDS"))

y2low <- as.data.frame.table(data2$y) %>%
  rename("individual" = "Var1",
         'interval' = "Var2",
         'site' = "Var3",
         "y" = "Freq")%>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual),
         interval = as.numeric(interval))

data3 <- readRDS(here("data_outputs",
                           "03_simulated",
                           "03_jags_input_data",
                           "Model3low_JAGS_data.RDS"))

y3low <- as.data.frame(data3$y) %>%
  rownames_to_column(var = "individual") %>%
  pivot_longer(2:last_col(),
               names_to = "site", 
               values_to = "y") %>%
  #mutate(site = str_sub(site, 4, nchar(site))) %>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual))

# GOF datasets ------------------------------------------------------------

gof1low <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod1low_GOFmodel.RDS'))

yrep1low <- as.data.frame.table(gof1low$sims.list$yrep) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "site" = "Var3",
         'yrep' = "Freq") %>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p1low <- as.data.frame.table(gof1low$sims.list$p) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "site" = "Var3",
         'p' = "Freq")%>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

gof2low <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod2low_GOFmodel.RDS'))

yrep2low <- as.data.frame.table(gof2low$sims.list$yrep) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         'interval' = "Var3",
         "site" = "Var4",
         'yrep' = "Freq")%>%
  mutate(site = as.numeric(site),
         interval = as.numeric(interval),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p2low <- as.data.frame.table(gof2low$sims.list$p.int) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         'interval' = "Var3",
         "site" = "Var4",
         'p' = "Freq")%>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual),
         interval = as.numeric(interval) ,
         iteration = as.numeric(iteration))

gof3low <- readRDS(here('data_outputs',
                             '03_simulated',
                             '04_posterior_summaries',
                             'mod3low_GOFmodel.RDS'))

yrep3.1low <- as.data.frame.table(gof3low$sims.list$yrep_1) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "site" = "Var3",
         'yrep' = "Freq")

yrep3.1lowid <- yrep3.1low %>%
  dplyr::select(iteration, individual, site)

yrep3.2low <- as.data.frame.table(gof3low$sims.list$yrep_2) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "site" = "Var3",
         'yrep' = "Freq") %>%
  anti_join(yrep3.1lowid, by = c("iteration", "individual", "site"))

yrep3low <- bind_rows(yrep3.1low, yrep3.2low)%>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p3.1low <- as.data.frame.table(gof3low$sims.list$pi1) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "site" = "Var3",
         'p' = "Freq")

p3.2low <- as.data.frame.table(gof3low$sims.list$pi2) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "site" = "Var3",
         'p' = "Freq") %>%
  anti_join(yrep3.1lowid, by = c("iteration", "individual", "site"))

p3low <- bind_rows(p3.1low, p3.2low)%>%
  mutate(site = as.numeric(site),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

# AUC ---------------------------------------------------------------------

iteration.num <- 1050

# Model 1 -----------------------------------------------------------------

resp1 <- as.vector(y1low$y)

mod1_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod1_AUC[i] <- AUC_JAGS2(p1low, 
                          iteration.num = i, 
                          resp = resp1)
}

mean1 <- as.data.frame(mod1_AUC) %>%
  summarise(mean = mean(mod1_AUC)) %>%
  as_vector()

(mod1_AUC_plot <- as.data.frame(mod1_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod1_AUC)) +
    geom_vline(xintercept = mean1, linetype = 2) +
    labs(title = "Total exposure") )

# Model 2 -----------------------------------------------------------------

resp2 <- as.vector(y2low$y)

mod2_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod2_AUC[i] <- AUC_JAGS2(p2low, 
                           iteration.num = i, 
                           resp = resp2)
}

mean2 <- as.data.frame(mod2_AUC) %>%
  summarise(mean = mean(mod2_AUC)) %>%
  as_vector()

(mod2_AUC_plot <- as.data.frame(mod2_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod2_AUC)) +
    geom_vline(xintercept = mean2, linetype = 2) +
    labs(title = "Interval") )

# Model 3 -----------------------------------------------------------------

resp3 <- as.vector(y3low$y)

mod3_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod3_AUC[i] <- AUC_JAGS2(p3low, 
                           iteration.num = i, 
                           resp = resp3)
}

mean3 <- as.data.frame(mod3_AUC) %>%
  summarise(mean = mean(mod3_AUC)) %>%
  as_vector()

(mod3_AUC_plot <- as.data.frame(mod3_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod3_AUC)) +
    geom_vline(xintercept = mean3, linetype = 2) +
    labs(title = "Custom") )

# Put together in one plot ------------------------------------------------

mod1AUC <- as.data.frame(mod1_AUC) %>%
  rename("AUC" = "mod1_AUC") %>%
  mutate(Model = "Total")

mod2AUC <- as.data.frame(mod2_AUC) %>%
  rename("AUC" = "mod2_AUC") %>%
  mutate(Model = "Interval")

mod3AUC <- as.data.frame(mod3_AUC) %>%
  rename("AUC" = "mod3_AUC") %>%
  mutate(Model = "Custom")

allAUC <- mod1AUC %>%
  bind_rows(mod2AUC, mod3AUC)

ggplot(allAUC) +
  geom_histogram(aes(x = AUC, fill = Model), 
                 alpha = 0.8) +
  geom_vline(xintercept = mean1, 
             linetype = 2,
             color = '#88419d') +
  geom_vline(xintercept = mean2, 
             linetype = 2,
             color = '#8c6bb1') +
  geom_vline(xintercept = mean3, 
             linetype = 2,
             color = '#8c96c6') +
  scale_fill_manual(values = c('#8c96c6',
                               '#8c6bb1',
                               '#88419d')) +
  labs(y = "Number of posterior samples")


# Export ------------------------------------------------------------------

saveRDS(allAUC, 
        here('data_outputs',
             '03_simulated',
             '05_for_plotting',
             'AUC_low.RDS'))


# Confusion matrices ------------------------------------------------------

conf_df1 <- yrep1low %>%
  left_join(y1low, by = c("individual", "site"))

conf_df2 <- yrep2low %>%
  left_join(y2low, by = c("individual", 
                          "site",
                          'interval'))

conf_df3 <- yrep3low %>%
  left_join(y3low, by = c("individual", "site"))

conf_summary_fun <- function(confusion_df){
  
  conf_summary <- confusion_df %>%
    mutate(category = case_when((yrep == 1 & y == 1) ~ "True Positive",
                                (yrep == 0 & y == 0) ~ "True Negative",
                                (yrep == 0 & y == 1) ~ "False Negative",
                                (yrep == 1 & y == 0) ~ "False Positive",
                                TRUE ~ NA_character_)) %>%
    filter(!is.na(category)) %>%
    group_by(iteration, category) %>%
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
  
  return(conf_summary)
}

conf_summary1 <- conf_summary_fun(conf_df1) %>%
  mutate(model = "Total")
conf_summary2 <- conf_summary_fun(conf_df2) %>%
  mutate(model = "Interval")
conf_summary3 <- conf_summary_fun(conf_df3) %>%
  mutate(model = "Custom") 

conf_all <- bind_rows(conf_summary1, conf_summary2,
                      conf_summary3)


ggplot(conf_all, aes(x = model, y = n, color = Classification)) +
  geom_boxplot() +
  facet_grid(Predicted ~ Actually, 
             scales = "free_x") +
  scale_y_log10() +
  scale_color_manual(values = c('#4d9221', 
                                '#c51b7d')) +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  labs(y= "Count")

saveRDS(conf_all, here("data_outputs", 
                           '01_whwonests',
                           '05_for_plotting',
                           'confusion_low.RDS'))







