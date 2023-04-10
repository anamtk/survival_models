# Model covariate results
# Ana Miller-ter Kuile
# March 24, 2022

# this script generates posterior distributions
# for the best-fitting survival model


# Load packages -----------------------------------------------------------


# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  'coda', #get posterior MCMC outputs
                  "bayesplot", #plot bayesian stuff
                  "tidybayes",#more bayesian plotting stuff
                  "mcmcplots", #posterior graphing
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

model_s <- readRDS(here("monsoon",
                        "9_29_22",
                        "survival",
                        "outputs",
                        "survival_mcmc_9_29.RDS"))


data_s <- read.csv(here("data_outputs",
                      "01_cleaning",
                      "03_nest_survival",
                      "Nest_survival_data.csv"))



data_s <- data_s %>%
  #some zero-day surveys that need to be deleted
  mutate(Julian_end = case_when(Julian_end == Julian_start ~ NA_integer_,
                                TRUE ~ Julian_end)) %>%
  filter(!is.na(Julian_end)) %>%
  #filter out multiple fledgling obs per nest
  filter(!(Stage == "F" & No_eggs == 999)) %>%
  group_by(Nest_ID) %>%
  #give a visit interval
  mutate(interval = 0:(n() - 1)) %>%
  mutate(prevStage = lag(Stage)) %>%
  ungroup() %>%
  #remove first survey
  filter(interval != 0) %>%  #1363
  filter(prevStage != "F") %>%
  #make incubating and laying one stage
  mutate(prevStage = case_when(prevStage == "E" ~ "Ex",
                               prevStage %in% c("I", "L") ~ "Eg",
                               prevStage == "N" ~ "Ne",
                               TRUE ~ NA_character_)) %>%
  filter(!prevStage == "Ex")


# Check p-values ----------------------------------------------------------


model_s$mean$zi
#nestheight (b[6])
#large trees b[9]
#precip b[13]
# forest cv b[14]
model_s$mean$zi.b1
#egg and nestling survival different
model_s$mean$zi.b2
model_s$mean$zi.b3
model_s$mean$zi.b4


# Medians and credible intervals for important variables ------------------
#nest heigth
model_s$q50$b[6]
model_s$q2.5$b[6]
model_s$q97.5$b[6]

#nest stage
model_s$q50$b1StageID
model_s$q2.5$b1StageID
model_s$q97.5$b1StageID

#large trees
model_s$q50$b[9]
model_s$q2.5$b[9]
model_s$q97.5$b[9]

#forest cv
model_s$q50$b[14]
model_s$q2.5$b[14]
model_s$q97.5$b[14]

#ppt
model_s$q50$b[13]
model_s$q2.5$b[13]
model_s$q97.5$b[13]

#overall intercept value
plogis(model_s$q50$b0)
plogis(model_s$q2.5$b0)
plogis(model_s$q97.5$b0)


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
