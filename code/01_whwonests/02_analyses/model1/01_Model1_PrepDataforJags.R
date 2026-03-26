# Nest survival data prep - full interval model1
# Ana Miller-ter Kuile
# March 31, 2023

# prepping data for the model of nest survival where data are
# derived only for the full survey interval

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

nests <- read.csv(here("data_outputs",
                       "01_whwonests",
                       "02_analysis_ready",
                       "full_survey_model_nest_data.csv"))


# Nest ID, number, and visit number ---------------------------------------
#get a vector of nest IDs
Nests <- unique(nests$Nest_ID)

# get the total count of nests
n.nests <- length(Nests)

# Random variables of forest, and year -------------------

#Year ID
Years <- nests %>%
  distinct(Year_located) %>%
  as_vector()

#number of years
n.years <- length(Years)

Forests <- nests %>% 
  distinct(Project_ID) %>%
  as_vector()

#number of forests
n.forests <- length(Forests)

# Create covariates -------------------------------------------------------

#the next sections create variables that encompass the random
# variables that group the data,
# and covariates that relate to the environment

# Random variables --------------------------------------------------------

#year as random effect - vector length of nests
Year <- nests %>%
  distinct(Nest_ID, Year_located) %>%
  dplyr::select(Year_located) %>%
  as_vector()

#make numeric for the model
Year.num <- nums(Year)

#forest random effect - vector of length of the 
#number of transects (hierarchically centered)
Forest <- nests %>%
  dplyr::select(Project_ID) %>%
  as_vector()

#make numeric for model
Forest.num <- nums(Forest)
        
# Nest and stand covariates -----------------------------------------------
#select all covariates on nest survival
nest_covs <- nests %>%
  dplyr::select(Nest_ID, Initiation_date,
                Nest_Ht, Orientation, 
                Trees_2550, Trees_50, 
                tmax, ppt) %>%
  #make a julian day
  mutate(Initiation_date = yday(as_date(Initiation_date))) %>%
  #set orientation to be 1 if N, -1 if south
  mutate(Orientation = cos(Orientation * (pi/180))) %>%
  mutate_if(is.numeric, scale)  #center and scale continous variables

#Nest-level covariates
NestHt <- as.vector(nest_covs$Nest_Ht)
InitDay <- as.vector(nest_covs$Initiation_date)
cosOrientation <- as.vector(nest_covs$Orientation) 
#Local covariates
Trees2550 <- as.vector(nest_covs$Trees_2550)
Trees50 <- as.vector(nest_covs$Trees_50)
#climate covarites
PPT <- as.vector(nest_covs$ppt)
Tmax <- as.vector(nest_covs$tmax)

# Exposure time -----------------------------------------------------------

#total amount of time each nest was surveyed
t <- as.vector(nests$t)

# Response vector ---------------------------------------------------------

y <- nests  %>%
  #set data to 1-0 for bernoulli
  mutate(survival = case_when(Fate_cat == "success" ~ 1,
                              Fate_cat == "failure" ~ 0,
                              TRUE ~ NA_real_)) %>%
  dplyr::select(survival) %>%
  #make a vector
  as_vector()

# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.nests = n.nests,
                 n.years = n.years,
                 n.forests = n.forests,
                 #Random effects IDs
                 Year.num = Year.num,
                 Forest.num = Forest.num,
                 #Nest-level covariates
                 NestHt = NestHt, 
                 cosOrientation = cosOrientation,
                 InitDay = InitDay, 
                 #Local-level covariates
                 Trees50 = Trees50,
                 Trees2550 = Trees2550, 
                #climate covariates
                 Tmax = Tmax,
                 PPT = PPT,
                 #dataset
                 y = y,
                 #exposure
                 t = t)

saveRDS(all_data, here("data_outputs", 
                       "01_whwonests",
                       '03_JAGS_input_data',
                      "mod1_JAGS_input_data.RDS"))

