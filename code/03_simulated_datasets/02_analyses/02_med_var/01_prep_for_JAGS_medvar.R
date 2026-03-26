# Prep simulated data for JAGS
# Ana Miller-ter Kuile
# May 17, 2023

# prep data from simulated datasets for jags


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  "jagsUI",
                  "coda",
                  "mcmcplots")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("code",
            "00_functions",
            "simulation_JAGS_data_list_funs.R"))

# Load data ---------------------------------------------------------------

#short daily data
med.day <- read.csv(here("data_outputs",
                         "03_simulated",
                         "02_analysis_ready",
                         "med_var_daily_data.csv"))

#short interval data
med.int <- read.csv(here("data_outputs",
                         "03_simulated",
                         "02_analysis_ready",
                         "med_var_interval_data.csv"))

#total data summarised
med.tot <- read.csv(here("data_outputs",
                         "03_simulated",
                         "02_analysis_ready",
                         "med_var_total_data.csv"))

# Model 1 -----------------------------------------------------------------

data1 <- exposure_model_dataprep_fun(med.tot)

saveRDS(data1, here("data_outputs",
                    "03_simulated",
                    "03_jags_input_data",
                    "Model1med_JAGS_data.RDS"))

# Model 2 -----------------------------------------------------------------

data2 <- interval_model_dataprep_fun(med.int)

#save it
saveRDS(data2, here("data_outputs",
                   "03_simulated",
                   "03_jags_input_data",
                   "Model2med_JAGS_data.RDS"))

# Model 3 -----------------------------------------------------------------

data3 <- custom_model_dataprep_fun(med.day)

saveRDS(data3, here("data_outputs",
                    "03_simulated",
                    "03_jags_input_data",
                    "Model3med_JAGS_data.RDS"))
