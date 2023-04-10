# Model covariate results
# Ana Miller-ter Kuile
# April 10, 2023

# this script generates posterior distributions
# for the survival model with interval-specific covariates
#but total end survival data


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

model3_sum <- readRDS(here("monsoon",
                           "01_whwonests",
                           "model3",
                           "outputs",
                           "model3_posterior_summary.RDS"))


# Posterior median and CI of all parameters -------------------------------

parms <- c("b[4]", "b[5]",
           "b[6]", "b[7]","b[8]",            
           "b[9]", "b[10]", "b[11]",          
           "b[12]", "b[13]", "b[14]",           
           "b[15]", "b[16]", "b[17]",  
           "b[18]",
           'b1StageID[2]',
           "b2TreatmentID[2]","b2TreatmentID[3]",
           "b2TreatmentID[4]", "b3SpeciesID[2]",  
           "b3SpeciesID[3]","b3SpeciesID[4]","b3SpeciesID[5]")

# "U", "B", "H", "HB"
#"PIPO", "Abies", "POTR5", "JUOC", "PSME"

mod3_est <- as.data.frame(model3_sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter %in% parms) %>%
  mutate(parameter = case_when(parameter == "b1StageID[2]" ~ 
                                 "Stage:Egg",
                               parameter == 'b2TreatmentID[2]' ~ 
                                 'TreatmentType:Burn',
                               parameter == "b2TreatmentID[3]" ~
                                 "TreatmentType:Harvest",
                               parameter == "b2TreatmentID[4]" ~
                                 "TreatmentType:Harvest+Burn",
                               parameter == "b3SpeciesID[2]" ~
                                 "NestTree:Abies",
                               parameter == "b3SpeciesID[3]" ~
                                 "NestTree:Aspen",
                               parameter == "b3SpeciesID[4]" ~
                                 "NestTree:Juniper",
                               parameter == "b3SpeciesID[5]" ~
                                 "NestTree:DougFir",
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
  mutate(Model = "Model3_CustomProbLogit")

write.csv(mod3_est, here("data_outputs",
                         "04_posterior_summaries",
                         "Model3_posteriors.csv"))


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