# Tree survival data prep - full interval model1
# Ana Miller-ter Kuile
# May 30, 2023

# prepping data for the model of tree survival where data are
# derived only for the full survey interval

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")

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
                      "total_tree_data.csv"))

# Plant ID, number --------------------------------------------------------

tree <- trees %>%
  distinct(CoreID)

#get a vector of plant IDS
Tree <- unique(trees$CoreID)

# get the total count of plants
n.trees <- length(Tree)


# Create covariates -------------------------------------------------------

#the next sections create variables that encompass the random
# variables that group the data,
# and covariates that relate to the environment

# Covariates --------------------------------------------------------------

covs <- trees %>%
  mutate(CoreID = as.factor(CoreID)) %>%
  mutate(Trt = as.factor(Trt)) %>%
  distinct(CoreID, Trt,
           BA, DBH, CanopyCover,
           meanVPD_Dry.Summer, meanVPD_Fall,
           meanVPD_Monsoon.Summer, meanVPD_Spring,
           meanVPD_Winter, 
           SWA_Dry.Summer, SWA_Fall,
           SWA_Monsoon.Summer, SWA_Spring,
           SWA_Winter) %>%
  mutate_if(is.numeric, scale) #center and scale continuous covariates

AvgDBH <- as.vector(covs$DBH) 
AvgBA <- as.vector(covs$BA)
CanopyCover <- as.vector(covs$CanopyCover)
VPD_ds <- as.vector(covs$meanVPD_Dry.Summer)
VPD_fa <- as.vector(covs$meanVPD_Fall)
VPD_ms <- as.vector(covs$meanVPD_Monsoon.Summer)
VPD_sp <- as.vector(covs$meanVPD_Spring)
VPD_wt <- as.vector(covs$meanVPD_Winter)
SWA_ds <- as.vector(covs$SWA_Dry.Summer)
SWA_fa <- as.vector(covs$SWA_Fall)
SWA_ms <- as.vector(covs$SWA_Monsoon.Summer)
SWA_sp <- as.vector(covs$SWA_Spring)
SWA_wt <- as.vector(covs$SWA_Winter)

covs %>%
  group_by(Trt) %>%
  tally()

TreatmentID <- nums(as.vector(covs$Trt))

n.trt <- length(unique(TreatmentID))

# Exposure time -----------------------------------------------------------

#total amount of time each plant was surveyed
t <- trees %>%
  rowwise() %>%
  mutate(t = end_year - start_year) %>%
  dplyr::select(t) %>%
  as_vector()

# Response vector ---------------------------------------------------------

y <- trees %>%
  mutate(status = case_when(response == "Live" ~ 1,
                            response == "Dead" ~ 0,
                            TRUE ~ NA_real_)) %>%
  dplyr::select(status) %>%
  #make a vector
  as_vector()

# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.trees = n.trees,
                 n.trt = n.trt,
                 #Covariates
                 AvgDBH = AvgDBH,
                 AvgBA = AvgBA,
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
                 #dataset
                 y = y,
                 #exposure
                 t = t)

saveRDS(all_data, here("data_outputs", 
                       "03_trees",
                       "03_JAGS_input_data",
                      "mod1_JAGS_input_data.RDS"))

saveRDS(all_data, 
        file = here("monsoon",
                    "03_trees",
                    "model1",
                    "inputs",
                    "mod1_JAGS_input_data.RDS"))

