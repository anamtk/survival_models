# Model covariate results
# Ana Miller-ter Kuile
# April 10, 2023

# this script generates posterior plots for all four models
#I've run so far


# Load packages -----------------------------------------------------------


# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "patchwork")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load posteriors ---------------------------------------------------------

m1 <- read.csv(here("data_outputs",
                    "04_posterior_summaries",
                    "Model1_posteriors.csv"))

m2 <- read.csv(here("data_outputs",
                    "04_posterior_summaries",
                    "Model2_posteriors.csv"))

m3 <- read.csv(here("data_outputs",
                    "04_posterior_summaries",
                    "Model3_posteriors.csv"))

# m3g <- read.csv(here("data_outputs",
#                      "04_posterior_summaries",
#                      "Model3gompit_posteriors.csv"))


# Merge data --------------------------------------------------------------

post <- bind_rows(m1, m2, m3)



# Plot --------------------------------------------------------------------

post <- post %>%
  mutate(parameter = factor(parameter,
                            levels = c("Stage:Egg", "Tmax", "Tmax^2",
                                       "PPT", "PPT^2",
                                       "NestHt", "NestOrientation",
                                       "InitDay", "NestTree:Abies",
                                       "NestTree:Aspen", "NestTree:Juniper",
                                       "NestTree:DougFir",
                                       "TreatmentType:Burn", 
                                       "TreatmentType:Harvest",
                                       "TreatmentType:Harvest+Burn",
                                       "SmTreeDens", "LgTreeDens",
                                       "PercPonderosa", "ForestCV",
                                       "Contagion", "NumOpenPatch",
                                       "PercHarvest", "PercBurn"))) %>%
  mutate(group = case_when(parameter %in% c("Stage:Egg", "Tmax", "Tmax^2",
                                            "PPT", "PPT^2") ~ "Interval",
                           parameter %in% c("NestHt", "NestOrientation",
                                            "InitDay", "NestTree:Abies",
                                            "NestTree:Aspen", "NestTree:Juniper",
                                            "NestTree:DougFir") ~ "Nest",
                           parameter %in% c("TreatmentType:Burn", 
                                            "TreatmentType:Harvest",
                                            "TreatmentType:Harvest+Burn",
                                            "PercHarvest", "PercBurn") ~ "Treatment",
                           parameter %in% c("SmTreeDens", "LgTreeDens",
                                            "PercPonderosa") ~ "Local",
                           parameter %in% c("ForestCV",
                                            "Contagion", "NumOpenPatch") ~ "Landscape")) %>%
  mutate(cat = case_when((X2.5. >0 & X50. > 0 & X97.5. > 0) ~ "sig",
                         (X2.5. < 0 & X50. < 0 & X97.5. < 0) ~ "sig",
                         TRUE ~ "notsig"))

post %>%
  #filter(Model != "Model3_CustomProbCloglog") %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(x = parameter, 
                      y= X50., 
                      ymin = X2.5.,
                      ymax = X97.5.,
                      shape = Model,
                      color = cat),
                  position = position_dodge2(width = 1)) +
  coord_flip() +
  scale_color_manual(values = c("black", "#af8dc3")) +
  facet_grid(group~., scales = "free") +
  theme_bw()
