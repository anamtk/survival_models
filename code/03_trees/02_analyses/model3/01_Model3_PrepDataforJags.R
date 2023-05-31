# Tree survival data prep
# November 1, 2021
# Ana Miller-ter Kuile

# prepping data for the model of tree survival

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

# Arrange by interval number for custom model -----------------------

trees2 <- trees %>%
  group_by(CoreID) %>%
  mutate(n.t = n()) %>%
  ungroup() %>%
  arrange(n.t)


# Count variables for indexing --------------------------------------------

# Numbers -----------------------------------------------------------------

#number with one interval only
n.trees1 <- trees2 %>%
  filter(n.t == 1) %>%
  tally() %>%
  as_vector()

#total number of plants
n.trees <- length(unique(trees2$CoreID))

#How many times did each nest get measured
#(number of intervals)
n.t <- as_vector(trees2$n.t)

n.trt <- length(unique(trees2$Trt))

# Covariates --------------------------------------------------------------

trees2 <- trees2 %>%  
  mutate(CoreID = factor(CoreID, levels = unique(CoreID))) %>%
  group_by(CoreID) %>%
  mutate(Visit_interval = 1:n()) %>%
  ungroup()


DBH <- int_cov2(variable = priorDBH)
BA <- int_cov2(variable = priorBA)
CanopyCover <- int_cov2(variable = CanopyCover)
VPD_ds <- int_cov2(variable = meanVPD_Dry.Summer)
VPD_fa <- int_cov2(variable = meanVPD_Fall)
VPD_ms <- int_cov2(variable = meanVPD_Monsoon.Summer)
VPD_sp <- int_cov2(variable = meanVPD_Spring)
VPD_wt <- int_cov2(variable = meanVPD_Winter)
SWA_ds <- int_cov2(variable = meanVPD_Dry.Summer)
SWA_fa <- int_cov2(variable = meanVPD_Fall)
SWA_ms <- int_cov2(variable = meanVPD_Monsoon.Summer)
SWA_sp <- int_cov2(variable = meanVPD_Spring)
SWA_wt <- int_cov2(variable = meanVPD_Winter)

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


# Response data vector ----------------------------------------------------

y <- trees2 %>%
  group_by(CoreID) %>%
  mutate(response = case_when(any(response == "Dead") ~ 0,
                              TRUE ~ 1)) %>%
  ungroup() %>%
  dplyr::select(response) %>%
  ungroup() %>%
  as_vector()

length(y[which(y == 1)])
length(y[which(y ==0)])

# Export as RDS -----------------------------------------------------------

all_data <- list(n.trees1 = n.trees1, 
                 n.trees = n.trees,
                 n.t = n.t,
                 n.trt = n.trt,
                 DBH = DBH,
                 BA = BA,
                 CanopyCover = CanopyCover,
                 VPD_ds = VPD_ds,
                 VPD_fa = VPD_fa,
                 VPD_ms = VPD_ms,
                 VPD_sp = VPD_sp,
                 VPD_wt = VPD_wt,
                 SWA_ds = SWA_ds,
                 SWA_fa = SWA_fa,
                 SWA_ms = SWA_ms,
                 SWA_sp = SWA_sp,
                 SWA_wt = SWA_wt,
                 TreatmentID = TreatmentID,
                 t = t,
                 y = y)


saveRDS(all_data, here("data_outputs", 
                       "03_trees",
                       '03_JAGS_input_data',
                       "mod3_JAGS_input_data.RDS"))

saveRDS(all_data, here("monsoon", 
                       "03_trees",
                       "model3",
                       "inputs",
                       "mod3_JAGS_input_data.RDS"))

