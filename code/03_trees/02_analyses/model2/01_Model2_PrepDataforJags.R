# Kelp survival data prep
# May 1, 2023
# Ana Miller-ter Kuile

# prepping data for the model of kelp surivival with interval data

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "lubridate", "glmmTMB")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("code",
            "00_functions",
            "tidy_functions.R"))

# Load data ---------------------------------------------------------------

trees <- read.csv(here("data_outputs",
                       "03_trees",
                       "02_analysis_ready",
                       "interval_tree_data.csv"))

# Variables of interest ---------------------------------------------------

# Count variables for indexing --------------------------------------------

n.trees <- length(unique(trees$CoreID))


#this should be a by-plant vector of number of visit intervals
n.t <- trees %>%
  mutate(CoreID = factor(CoreID, levels = unique(CoreID))) %>%
  group_by(CoreID) %>%
  summarise(n.t = n()) %>%
  dplyr::select(n.t) %>%
  as_vector()

n.trt <- length(unique(trees$Trt))



# Covariates --------------------------------------------------------------

trees2 <- trees %>%  
  mutate(CoreID = factor(CoreID, levels = unique(CoreID))) %>%
  group_by(CoreID) %>%
  mutate(Visit_interval = 1:n()) %>%
  ungroup()
  

DBH <- int_cov2(variable = priorDBH)
BA <- int_cov2(variable = priorBA)
CanopyCover <- int_cov2(variable = CanopyCover)
maxVPD <- int_cov2(variable = maxVPD)
minSWA <- int_cov2(variable = minSWA)

TreatmentID <- trees2 %>%
  dplyr::select(CoreID, Visit_interval, Trt) %>%
  mutate(Trt = as.numeric(as.factor(Trt))) %>%
  pivot_wider(names_from = "Visit_interval",
              values_from = "Trt")  %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()
  
# Time interval matrix ----------------------------------------------------

# this should be a plant x interval matrix
t <- trees2 %>%
  rowwise() %>%
  mutate(t = Year - priorYear) %>%
  dplyr::select(CoreID, Visit_interval, t) %>%
  pivot_wider(names_from = Visit_interval,
              values_from = t) %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()


# Response data -----------------------------------------------------------

y <- trees2 %>%
  mutate(response = case_when(response == "Live" ~ 1,
                              response == "Dead" ~ 0)) %>%
  dplyr::select(CoreID, Visit_interval, response) %>%
  pivot_wider(names_from = Visit_interval,
              values_from = response) %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()


# Compile and export ------------------------------------------------------

all_data <- list(n.trees = n.trees,
                 n.t = n.t,
                 n.trt = n.trt,
                 DBH = DBH,
                 BA = BA,
                 CanopyCover = CanopyCover,
                 maxVPD = maxVPD,
                 minSWA = minSWA,
                 TreatmentID = TreatmentID,
                 t = t,
                 y = y)

saveRDS(all_data, here("data_outputs", 
                       "03_trees",
                       "03_JAGS_input_data",
                       "mod2_JAGS_input_data.RDS"))

saveRDS(all_data, 
        file = here("monsoon",
                    "03_trees",
                    "model2",
                    "inputs",
                    "mod2_JAGS_input_data.RDS"))





