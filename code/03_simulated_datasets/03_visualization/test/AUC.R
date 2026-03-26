

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

data1_1 <- readRDS(here("data_outputs",
                      "03_simulated",
                      "03_jags_input_data",
                      'test',
                      "Model1low_JAGS_data_1.RDS"))


y1_1 <- as.data.frame(data1_1$y)%>%
  rownames_to_column(var = "individual") %>%
  pivot_longer(2:last_col(),
               names_to = "dataset", 
               values_to = "y") %>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual))

data1_2 <- readRDS(here("data_outputs",
                        "03_simulated",
                        "03_jags_input_data",
                        'test',
                        "Model1low_JAGS_data_2.RDS"))

y1_2 <- as.data.frame(data1_2$y)%>%
  rownames_to_column(var = "individual") %>%
  pivot_longer(2:last_col(),
               names_to = "dataset", 
               values_to = "y") %>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual))


data2_1 <- readRDS(here("data_outputs",
                           "03_simulated",
                           "03_jags_input_data",
                        'test',
                        "Model2low_JAGS_data_1.RDS"))

y2_1 <- as.data.frame.table(data2_1$y) %>%
  rename("individual" = "Var1",
         'interval' = "Var2",
         'dataset' = "Var3",
         "y" = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         interval = as.numeric(interval))

data2_2 <- readRDS(here("data_outputs",
                        "03_simulated",
                        "03_jags_input_data",
                        'test',
                        "Model2low_JAGS_data_2.RDS"))

y2_2 <- as.data.frame.table(data2_2$y) %>%
  rename("individual" = "Var1",
         'interval' = "Var2",
         'dataset' = "Var3",
         "y" = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         interval = as.numeric(interval))

data3_1 <- readRDS(here("data_outputs",
                           "03_simulated",
                           "03_jags_input_data",
                        'test',
                           "Model3low_JAGS_data_1.RDS"))


y3_1 <- as.data.frame(data3_1$y) %>%
  rownames_to_column(var = "individual") %>%
  pivot_longer(2:last_col(),
               names_to = "dataset", 
               values_to = "y") %>%
  mutate(dataset = str_sub(dataset, 4, nchar(dataset))) %>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual))

data3_2 <- readRDS(here("data_outputs",
                        "03_simulated",
                        "03_jags_input_data",
                        'test',
                        "Model3low_JAGS_data_2.RDS"))


y3_2 <- as.data.frame(data3_2$y) %>%
  rownames_to_column(var = "individual") %>%
  pivot_longer(2:last_col(),
               names_to = "dataset", 
               values_to = "y") %>%
  mutate(dataset = str_sub(dataset, 4, nchar(dataset))) %>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual))


# GOF datasets ------------------------------------------------------------

gof1_1 <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                       'test',
                        'mod1low1_GOFmodel.RDS'))

yrep1_1 <- as.data.frame.table(gof1_1$sims.list$yrep) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'yrep' = "Freq") %>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p1_1 <- as.data.frame.table(gof1_1$sims.list$p) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'p' = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

#dataset2 - model 1
gof1_2 <- readRDS(here('data_outputs',
                       '03_simulated',
                       '04_posterior_summaries',
                       'test',
                       'mod1low2_GOFmodel.RDS'))

yrep1_2 <- as.data.frame.table(gof1_2$sims.list$yrep) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'yrep' = "Freq") %>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p1_2 <- as.data.frame.table(gof1_2$sims.list$p) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'p' = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

#model 2 dataset 1

gof2_1 <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                       'test',
                        'mod2low1_GOFmodel.RDS'))

yrep2_1 <- as.data.frame.table(gof2_1$sims.list$yrep) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         'interval' = "Var3",
         "dataset" = "Var4",
         'yrep' = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         interval = as.numeric(interval),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p2_1 <- as.data.frame.table(gof2_1$sims.list$p.int) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         'interval' = "Var3",
         "dataset" = "Var4",
         'p' = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         interval = as.numeric(interval) ,
         iteration = as.numeric(iteration))

#model 2 dataset 2
gof2_2 <- readRDS(here('data_outputs',
                       '03_simulated',
                       '04_posterior_summaries',
                       'test',
                       'mod2low2_GOFmodel.RDS'))

yrep2_2 <- as.data.frame.table(gof2_2$sims.list$yrep) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         'interval' = "Var3",
         "dataset" = "Var4",
         'yrep' = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         interval = as.numeric(interval),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p2_2 <- as.data.frame.table(gof2_2$sims.list$p.int) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         'interval' = "Var3",
         "dataset" = "Var4",
         'p' = "Freq")%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         interval = as.numeric(interval) ,
         iteration = as.numeric(iteration))

#model 3 dataset 1

gof3_1 <- readRDS(here('data_outputs',
                             '03_simulated',
                             '04_posterior_summaries',
                       'test',
                             'mod3low2_GOFmodel.RDS'))

yrep3_1.1 <- as.data.frame.table(gof3_1$sims.list$yrep_1) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'yrep' = "Freq")

#get IDs to anti join with yrep 2
yrep3_1_ids <- yrep3_1.1 %>%
  dplyr::select(iteration, individual, dataset)

yrep3_1.2 <-  as.data.frame.table(gof3_1$sims.list$yrep_2) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'yrep' = "Freq") %>%
  anti_join(yrep3_1_ids, 
            by = c("iteration", "individual", 'dataset'))

yrep3_1 <- bind_rows(yrep3_1.1, yrep3_1.2)%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p3_1.1 <- as.data.frame.table(gof3_1$sims.list$pi1) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'p' = "Freq")

p3_1.2 <- as.data.frame.table(gof3_1$sims.list$pi2) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'p' = "Freq") %>%
  anti_join(yrep3_1_ids, 
            by = c("iteration", "individual", 'dataset'))

p3_1 <- bind_rows(p3_1.1, p3_1.2)%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

#model 3 dataset 2

gof3_2 <- readRDS(here('data_outputs',
                       '03_simulated',
                       '04_posterior_summaries',
                       'test',
                       'mod3low2_GOFmodel.RDS'))

yrep3_2.1 <- as.data.frame.table(gof3_2$sims.list$yrep_1) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'yrep' = "Freq")

yrep3_2_ids <- yrep3_2.1 %>%
  dplyr::select(iteration, individual, dataset)

yrep3_2.2 <- as.data.frame.table(gof3_2$sims.list$yrep_2) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'yrep' = "Freq")%>%
  anti_join(yrep3_2_ids, 
            by = c("iteration", "individual", 'dataset'))

yrep3_2 <- bind_rows(yrep3_2.1, yrep3_2.2)%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

p3_2.1 <- as.data.frame.table(gof3_2$sims.list$pi1) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'p' = "Freq")

p3_2.2 <- as.data.frame.table(gof3_2$sims.list$pi2) %>%
  rename('iteration' = "Var1",
         "individual" = "Var2",
         "dataset" = "Var3",
         'p' = "Freq")%>%
  anti_join(yrep3_2_ids, 
            by = c("iteration", "individual", 'dataset'))

p3_2 <- bind_rows(p3_2.1, p3_2.2)%>%
  mutate(dataset = as.numeric(dataset),
         individual = as.numeric(individual),
         iteration = as.numeric(iteration))

# AUC ---------------------------------------------------------------------

iteration.num <- 1050

# Model 1 -----------------------------------------------------------------

#dataset 1
resp1_1 <- as.vector(y1_1$y)

mod1_1_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod1_1_AUC[i] <- AUC_JAGS2(p1_1, 
                          iteration.num = i, 
                          resp = resp1_1)
}

mean1_1 <- as.data.frame(mod1_1_AUC) %>%
  summarise(mean = mean(mod1_1_AUC)) %>%
  as_vector()

(mod1_1_AUC_plot <- as.data.frame(mod1_1_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod1_1_AUC)) +
    geom_vline(xintercept = mean1_1, linetype = 2) +
    labs(title = "Total exposure") )

#dataset 2
resp1_2 <- as.vector(y1_2$y)

mod1_2_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod1_2_AUC[i] <- AUC_JAGS2(p1_2, 
                             iteration.num = i, 
                             resp = resp1_2)
}

mean1_2 <- as.data.frame(mod1_2_AUC) %>%
  summarise(mean = mean(mod1_2_AUC)) %>%
  as_vector()

(mod1_2_AUC_plot <- as.data.frame(mod1_2_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod1_2_AUC)) +
    geom_vline(xintercept = mean1_2, linetype = 2) +
    labs(title = "Total exposure") )


# Model 2 -----------------------------------------------------------------

#dataset 1
resp2_1 <- as.vector(y2_1$y)

mod2_1_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod2_1_AUC[i] <- AUC_JAGS2(p2_1, 
                           iteration.num = i, 
                           resp = resp2_1)
}

mean2_1 <- as.data.frame(mod2_1_AUC) %>%
  summarise(mean = mean(mod2_1_AUC)) %>%
  as_vector()

(mod2_1_AUC_plot <- as.data.frame(mod2_1_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod2_1_AUC)) +
    geom_vline(xintercept = mean2_1, linetype = 2) +
    labs(title = "Interval") )


#dataset2 
resp2_2 <- as.vector(y2_2$y)

mod2_2_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod2_2_AUC[i] <- AUC_JAGS2(p2_2, 
                                iteration.num = i, 
                                resp = resp2_2)
}

mean2_2 <- as.data.frame(mod2_2_AUC) %>%
  summarise(mean = mean(mod2_2_AUC)) %>%
  as_vector()

(mod2_2_AUC_plot <- as.data.frame(mod2_2_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod2_2_AUC)) +
    geom_vline(xintercept = mean2_2, linetype = 2) +
    labs(title = "Interval") )


# Model 3 -----------------------------------------------------------------

#dataset1 
resp3_1 <- as.vector(y3_1$y)

mod3_1_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod3_1_AUC[i] <- AUC_JAGS2(p3_1, 
                           iteration.num = i, 
                           resp = resp3_1)
}

mean3_1 <- as.data.frame(mod3_1_AUC) %>%
  summarise(mean = mean(mod3_1_AUC)) %>%
  as_vector()

(mod3_1_AUC_plot <- as.data.frame(mod3_1_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod3_1_AUC)) +
    geom_vline(xintercept = mean3_1, linetype = 2) +
    labs(title = "Custom") )

#dataset 2
resp3_2 <- as.vector(y3_2$y)

mod3_2_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod3_2_AUC[i] <- AUC_JAGS2(p3_2, 
                                iteration.num = i, 
                                resp = resp3_2)
}

mean3_2 <- as.data.frame(mod3_2_AUC) %>%
  summarise(mean = mean(mod3_2_AUC)) %>%
  as_vector()

(mod3_2_AUC_plot <- as.data.frame(mod3_2_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod3_2_AUC)) +
    geom_vline(xintercept = mean3_2, linetype = 2) +
    labs(title = "Custom") )

# Put together in one plot ------------------------------------------------

mod1_1AUC <- as.data.frame(mod1_1_AUC) %>%
  rename("AUC" = "mod1_1_AUC") %>%
  mutate(Model = "Total") %>%
  mutate(Dataset = "1")

mod1_2AUC <- as.data.frame(mod1_2_AUC) %>%
  rename("AUC" = "mod1_2_AUC") %>%
  mutate(Model = "Total") %>%
  mutate(Dataset = "2")

mod2_1AUC <- as.data.frame(mod2_1_AUC) %>%
  rename("AUC" = "mod2_1_AUC") %>%
  mutate(Model = "Interval")%>%
  mutate(Dataset = "1")

mod2_2AUC <- as.data.frame(mod2_2_AUC) %>%
  rename("AUC" = "mod2_2_AUC") %>%
  mutate(Model = "Interval")%>%
  mutate(Dataset = "2")

mod3_1AUC <- as.data.frame(mod3_1_AUC) %>%
  rename("AUC" = "mod3_1_AUC") %>%
  mutate(Model = "Custom") %>%
  mutate(Dataset = "1")

mod3_2AUC <- as.data.frame(mod3_2_AUC) %>%
  rename("AUC" = "mod3_2_AUC") %>%
  mutate(Model = "Custom") %>%
  mutate(Dataset = "2")

allAUC <- mod1_1AUC %>%
  bind_rows(mod1_2AUC, mod2_1AUC, 
            mod2_2AUC, mod3_1AUC,
            mod3_2AUC)

ggplot(allAUC) +
  geom_histogram(aes(x = AUC, fill = Model), 
                 alpha = 0.8) +
  # geom_vline(xintercept = mean1, 
  #            linetype = 2,
  #            color = '#88419d') +
  # geom_vline(xintercept = mean2, 
  #            linetype = 2,
  #            color = '#8c6bb1') +
  # geom_vline(xintercept = mean3, 
  #            linetype = 2,
  #            color = '#8c96c6') +
  scale_fill_manual(values = c('#8c96c6',
                               '#8c6bb1',
                               '#88419d')) +
  labs(y = "Number of posterior samples") +
  facet_wrap(~Dataset)


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

saveRDS(conf_summary, here("data_outputs", 
                           '01_whwonests',
                           '05_for_plotting',
                           'confusion_mod1.RDS'))







