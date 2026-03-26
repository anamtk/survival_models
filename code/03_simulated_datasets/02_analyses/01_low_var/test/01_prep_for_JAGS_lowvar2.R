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

#daily data
low.day <- read.csv(here("data_outputs",
                         "03_simulated",
                         "02_analysis_ready",
                         'test',
                         "low_var_daily_data_2.csv"))

#interval data
low.int <- read.csv(here("data_outputs",
                         "03_simulated",
                         "02_analysis_ready",
                         'test',
                         "low_var_interval_data_2.csv"))

#total data summarised
low.tot <- read.csv(here("data_outputs",
                         "03_simulated",
                         "02_analysis_ready",
                         'test',
                         "low_var_total_data_2.csv"))

# Model 1 -----------------------------------------------------------------

data1 <- exposure_model_dataprep_fun(low.tot)

saveRDS(data1, here("data_outputs",
                    "03_simulated",
                    "03_jags_input_data",
                    'test',
                    "Model1low_JAGS_data_2.RDS"))

# Model 2 -----------------------------------------------------------------

data2 <- interval_model_dataprep_fun(low.int)

#save it
saveRDS(data2, here("data_outputs",
                   "03_simulated",
                   "03_jags_input_data",
                   'test',
                   "Model2low_JAGS_data_2.RDS"))

# Model 3 -----------------------------------------------------------------

data3 <- custom_model_dataprep_fun(low.day)

saveRDS(data3, here("data_outputs",
                    "03_simulated",
                    "03_jags_input_data",
                    'test',
                    "Model3low_JAGS_data_2.RDS"))



