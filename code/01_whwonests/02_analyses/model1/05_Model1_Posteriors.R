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

parms <- c("b[3]", "b[4]", "b[5]",
           "b[6]", "b[7]","b[8]",            
           "b[9]", "b[10]", "b[11]",          
           "b[12]", "b[13]", "b[14]",           
           "b[15]", "b[16]", "b[17]",           
           "b1TreatmentID[2]","b1TreatmentID[3]",
           "b1TreatmentID[4]", "b2SpeciesID[2]",  
           "b2SpeciesID[3]","b2SpeciesID[4]","b2SpeciesID[5]")

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
                               parameter == "b[3]" ~ "NestHt",
                               parameter == "b[4]" ~ "NestOrientation",
                               parameter == "b[5]" ~ "InitDay",
                               parameter == "b[6]" ~ "LgTreeDens",
                               parameter == "b[7]" ~ "SmTreeDens",
                               parameter == "b[8]" ~ "PercPonderosa",
                               parameter == "b[9]" ~ "Tmax",
                               parameter == "b[10]" ~ "Tmax^2",
                               parameter == "b[11]" ~ "PPT",
                               parameter == "b[12]" ~ "PPT^2",
                               parameter == "b[13]" ~ "ForestCV",
                               parameter == "b[14]" ~ "Contagion",
                               parameter == "b[15]" ~ "NumOpenPatch",
                               parameter == "b[16]" ~ "PercHarvest",
                               parameter == "b[17]" ~ "PercBurn",
                               TRUE ~ parameter)) %>%
  mutate(Model = "Model1_TotalExposure")
  

write.csv(mod1_est, here("data_outputs",
               "04_posterior_summaries",
               "Model1_posteriors.csv"))
# Sources of variation ----------------------------------------------------

model_s$q50$sig.nest
model_s$q2.5$sig.nest
model_s$q97.5$sig.nest

model_s$q50$sig.transect
model_s$q2.5$sig.transect
model_s$q97.5$sig.transect

model_s$q50$sig.year
model_s$q2.5$sig.year
model_s$q97.5$sig.year
