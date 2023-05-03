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

model3_sum <- readRDS(here("monsoon",
                           "02_kelp",
                           "model3",
                           "outputs",
                           "model3_posterior_summary.RDS"))


# Posterior median and CI of all parameters -------------------------------

parms <- c("b[1]", "b[2]",
           "b[3]", "b[4]","b[5]",  
           'b6[2]', 'b6[3]',
           'b6[4]', 'b6[5]')

# "B", "BO",
# "C", "S", "SS"

mod3_est <- as.data.frame(model3_sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter %in% parms) %>%
  mutate(parameter = case_when(parameter == "b[1]" ~ "Sea surface temp",
                               parameter == "b[2]" ~ "Wave power",
                               parameter == "b[3]" ~ "Stipe number",
                               parameter == "b[4]" ~ "Holdfast diameter",
                               parameter == "b[5]" ~ "Plant depth",
                               parameter == 'b6[2]' ~ "Substrate = boulder",
                               parameter == "b6[3]" ~ "Substrate = cobble",
                               parameter == "b6[4]" ~ "Substrate = sand",
                               parameter == "b6[5]" ~ "Substrate = shallow sand",
                               TRUE ~ parameter)) %>%
  mutate(Model = "Model3_CustomProb")

write.csv(mod3_est, here("data_outputs",
                         "02_kelp",
                         "04_posterior_summaries",
                         "Model3_posteriors.csv"))


#