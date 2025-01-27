#Look at predictive accruacy of each model
#Ana Miller-ter Kuile
#June 12, 2023  

#this script pulls out the posterior ranges of each parameter to see 
# how often the model predicts a value in the range of the correct value

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("code",
            "00_functions",
            "simulation_functions.R"))

theme_set(theme_bw())

# Low var -----------------------------------------------------------------

# Load datasets -----------------------------------------------------------

mod1low <- readRDS(here("monsoon",
                        "simulated",
                        "lowvar",
                        "outputs",
                        'mod1_lowvar_sum.RDS'))

mod2low <- readRDS(here("monsoon",
                        "simulated",
                        "lowvar",
                        "outputs",
                        'mod2_lowvar_sum.RDS'))

mod3low <- readRDS(here("monsoon",
                        "simulated",
                        "lowvar",
                        "outputs",
                        'mod3_lowvar_sum.RDS'))

# Med var -----------------------------------------------------------------

# Load datasets -----------------------------------------------------------

mod1med <- readRDS(here("monsoon",
                        "simulated",
                        "medvar",
                        "outputs",
                        'mod1_medvar_sum.RDS'))

mod2med <- readRDS(here("monsoon",
                        "simulated",
                        "medvar",
                        "outputs",
                        'mod2_medvar_sum.RDS'))

mod3med <- readRDS(here("monsoon",
                        "simulated",
                        "medvar",
                        "outputs",
                        'mod3_medvar_sum.RDS'))

# High Var ----------------------------------------------------------------

# Load datasets -----------------------------------------------------------

mod1high <- readRDS(here("monsoon",
                        "simulated",
                        "highvar",
                        "outputs",
                        'mod1_highvar_sum.RDS'))

mod2high <- readRDS(here("monsoon",
                        "simulated",
                        "highvar",
                        "outputs",
                        'mod2_highvar_sum.RDS'))

mod3high <- readRDS(here("monsoon",
                        "simulated",
                        "highvar",
                        "outputs",
                        'mod3_highvar_sum.RDS'))

# Prep data ---------------------------------------------------------------

#b0 <- 0.5 #gets mean value survival fairly high
#b1 <- 0.5 #this value makes sure that ps is always >0.5 in the range of x

post_perc_fun(mod1 = mod1low,
              mod2 = mod2low,
              mod3 = mod3low,
              parm = 0.5,
              beta = "b1")

post_perc_fun(mod1 = mod1med,
              mod2 = mod2med,
              mod3 = mod3med,
              parm = 0.5,
              beta = "b1")

post_perc_fun(mod1 = mod1high,
              mod2 = mod2high,
              mod3 = mod3high,
              parm = 0.5,
              beta = "b1")


post_perc_fun(mod1 = mod1low,
              mod2 = mod2low,
              mod3 = mod3low,
              parm = 0.5,
              beta = "b0")

post_perc_fun(mod1 = mod1med,
              mod2 = mod2med,
              mod3 = mod3med,
              parm = 0.5,
              beta = "b0")

post_perc_fun(mod1 = mod1high,
              mod2 = mod2high,
              mod3 = mod3high,
              parm = 0.5,
              beta = "b0")

