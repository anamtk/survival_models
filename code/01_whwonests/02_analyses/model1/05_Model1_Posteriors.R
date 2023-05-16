# Model covariate results
# Ana Miller-ter Kuile
# April 10, 2023

# this script generates posterior distributions
# for the survival model with full interval


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

source(here::here("code", 
                  "00_functions",
                  "plot_functions.R"))

source(here::here("code",
                  "00_functions",
                  "tidy_functions.R"))

# Load data ---------------------------------------------------------------

model1_sum <- readRDS(here("monsoon",
                           "01_whwonests",
                           "model1",
                           "outputs",
                           "model1_posterior_summary.RDS"))


# Posterior median and CI of all parameters -------------------------------

parms <- c("b[4]", "b[5]",
           "b[6]", "b[7]","b[8]",            
           "b[9]", "b[10]", "b[11]",          
           "b[12]", "b[13]", "b[14]",           
           "b[15]", "b[16]", "b[17]",  
           'b[18]',
           "b1TreatmentID[2]","b1TreatmentID[3]",
           "b1TreatmentID[4]", "b2SpeciesID[2]",  
           "b2SpeciesID[3]","b2SpeciesID[4]","b2SpeciesID[5]",
           'b3StageID[2]')

# "U", "B", "H", "HB"
#"PIPO", "Abies", "POTR5", "JUOC", "PSME"

mod1_est <- as.data.frame(model1_sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter %in% parms) %>%
  mutate(parameter = case_when(parameter == 'b1TreatmentID[2]' ~ 
                                 'TreatmentType:Burn',
                               parameter == "b1TreatmentID[3]" ~
                                 "TreatmentType:Harvest",
                               parameter == "b1TreatmentID[4]" ~
                                 "TreatmentType:Harvest+Burn",
                               parameter == "b2SpeciesID[2]" ~
                                 "NestTree:Abies",
                               parameter == "b2SpeciesID[3]" ~
                                 "NestTree:Aspen",
                               parameter == "b2SpeciesID[4]" ~
                                 "NestTree:Juniper",
                               parameter == "b2SpeciesID[5]" ~
                                 "NestTree:DougFir",
                               parameter == "b3StageID[2]" ~
                                 "Stage:Egg",
                               parameter == "b[4]" ~ "NestHt",
                               parameter == "b[5]" ~ "NestOrientation",
                               parameter == "b[6]" ~ "InitDay",
                               parameter == "b[7]" ~ "LgTreeDens",
                               parameter == "b[8]" ~ "SmTreeDens",
                               parameter == "b[9]" ~ "PercPonderosa",
                               parameter == "b[10]" ~ "Tmax",
                               parameter == "b[11]" ~ "Tmax^2",
                               parameter == "b[12]" ~ "PPT",
                               parameter == "b[13]" ~ "PPT^2",
                               parameter == "b[14]" ~ "ForestCV",
                               parameter == "b[15]" ~ "Contagion",
                               parameter == "b[16]" ~ "NumOpenPatch",
                               parameter == "b[17]" ~ "PercHarvest",
                               parameter == "b[18]" ~ "PercBurn",
                               TRUE ~ parameter)) %>%
  mutate(Model = "Model1_TotalExposure")
  

write.csv(mod1_est, here("data_outputs",
                         "01_whwonests",
               "04_posterior_summaries",
               "Model1_posteriors.csv"))


# P-values ----------------------------------------------------------------

p_values1 <- as.data.frame(model1_sum$statistics) %>%
  dplyr::select(Mean) %>%
  rownames_to_column(var = "zvalue") %>%
  filter(str_detect(zvalue, "z")) %>%
  mutate(p = case_when(Mean >= 0.5 ~ (1-Mean), #I Think these are 1-tailed p-values but check with Kiona
                       Mean < 0.5 ~ (1 - (1-Mean)))) %>%
  mutate(p_cat = case_when(p <= 0.05 ~ "s",
                           TRUE ~ "ns")) %>%
  mutate(direction = case_when(Mean >= 0.5 ~ 'positive', 
                       Mean < 0.5 ~ "negative")) %>%
  mutate(p_dir = case_when((p_cat == "s" & direction == "positive") ~ p,
                           (p_cat == "s" & direction == "negative") ~ -p,
                           TRUE ~ NA_real_)) %>%
  filter(!zvalue %in% c("z.b1[1]", "z.b2[1]", 'z.b3[1]')) %>%
  mutate(zvalue = case_when(zvalue == 'z.b1[2]' ~ 
                                 'TreatmentType:Burn',
                               zvalue == "z.b1[3]" ~
                                 "TreatmentType:Harvest",
                               zvalue == "z.b1[4]" ~
                                 "TreatmentType:Harvest+Burn",
                               zvalue == "z.b2[2]" ~
                                 "NestTree:Abies",
                               zvalue == "z.b2[3]" ~
                                 "NestTree:Aspen",
                               zvalue == "z.b2[4]" ~
                                 "NestTree:Juniper",
                               zvalue == "z.b2[5]" ~
                                 "NestTree:DougFir",
                               zvalue == "z.b3[2]" ~
                                 "Stage:Egg",
                               zvalue == "z[4]" ~ "NestHt",
                               zvalue == "z[5]" ~ "NestOrientation",
                               zvalue == "z[6]" ~ "InitDay",
                               zvalue == "z[7]" ~ "LgTreeDens",
                               zvalue == "z[8]" ~ "SmTreeDens",
                               zvalue == "z[9]" ~ "PercPonderosa",
                               zvalue == "z[10]" ~ "Tmax",
                               zvalue == "z[11]" ~ "Tmax^2",
                               zvalue == "z[12]" ~ "PPT",
                               zvalue == "z[13]" ~ "PPT^2",
                               zvalue == "z[14]" ~ "ForestCV",
                               zvalue == "z[15]" ~ "Contagion",
                               zvalue == "z[16]" ~ "NumOpenPatch",
                               zvalue == "z[17]" ~ "PercHarvest",
                               zvalue == "z[18]" ~ "PercBurn",
                               TRUE ~ zvalue)) %>%
  mutate(Model = "Model1_TotalExposure")

p_values1 %>%
  filter(!is.na(p_dir)) %>%
  ggplot(aes(x = reorder(zvalue, p), y = p_dir)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  coord_flip()


write.csv(p_values1, here("data_outputs",
                         "01_whwonests",
                         "04_posterior_summaries",
                         "Model1_pvalues.csv"))
