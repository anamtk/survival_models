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

Nests <- nests %>%
  distinct(Nest_ID)

#get a vector of nest IDs
Nests <- unique(nests$Nest_ID)

# get the total count of nests
n.nests <- length(Nests)

# Random variables of transect, forest, and year -------------------

#Get the names of all the transects
Transects <- nests %>%
  distinct(Transect_ID2) %>%
  as_vector()

#number of transects
n.transects <- length(Transects)

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

#transect random effect (vector length of number of nest points)
Transect <- nests %>%
  distinct(Nest_ID, 
           Transect_ID2,
           Project_ID) %>%
  dplyr::select(Transect_ID2) %>%
  as_vector()

#make numeric for the model
Transect.num <- nums(Transect)

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
  distinct(Transect_ID2, 
           Project_ID) %>%
  dplyr::select(Project_ID) %>%
  as_vector()

#make numeric for model
Forest.num <- nums(Forest)

# MIssing data variables --------------------------------------------------

#vector length of number of nests of forest ID for
#imputing the climate data

Forest1 <- nests %>%
  distinct(Nest_ID,
           Transect_ID2,
           Project_ID) %>%
  dplyr::select(Project_ID) %>%
  as_vector()

#get taht in numeric for model
Forest.ID <- nums(Forest1)
           
# Nest and stand covariates -----------------------------------------------
#select all covariates on nest survival
nest_covs <- nests %>%
  dplyr::select(Nest_ID, Nest_Ht, 
           Tree_sp, cosOrientation, Init_day, Trt_cat, prevStage,
           Trees_2550, Trees_50, pPIPO,
           a1000_areacv2, a1000_contag, a1000_np1, a1000_Ha, a1000_RxBu,
           meanTmax_C, meanPpt_mm) %>% 
  mutate_if(is.numeric, scale)  #center and scale continous variables

#categorical covariates
#set the base level to be the one with the 
#most observations

#Treatment covariate
TreatmentID <- nest_covs %>%
  mutate(Trt_cat = factor(Trt_cat, levels = c("U", "B", "H", "HB"))) %>%
  dplyr::select(Trt_cat) %>%
  as_vector() %>%
  nums()

n.trt <- length(unique(as.factor(TreatmentID)))
#levels:
#1 = untreated
#2 = burn
#3 = harvest
#4 = harvest burn

## Species categorical effect
SpeciesID <- nest_covs %>%
  mutate(Tree_sp = factor(Tree_sp, levels = c("PIPO", "Abies", "POTR5",
                                              "JUOC", "PSME"))) %>%
  dplyr::select(Tree_sp) %>%
  as_vector() %>%
  nums()
  
n.species <- length(unique(as.factor(SpeciesID)))

#levels(as.factor(nest_covs$Tree_sp))
#1 = PIPO
#2 = Abies
#3 = POTR5
#4 = JUOC
#5 = PSME

#stage categorical effect
nests %>%
  group_by(prevStage) %>%
  tally()

StageID <- nests %>%
  mutate(StageID = case_when(prevStage == "N" ~ 1,
                             prevStage %in% c("I", "L") ~ 2,
                             TRUE ~ NA_real_)) %>%
  dplyr::select(StageID) %>%
  as_vector() %>%
  nums()

n.stages <- 2

#Nest-level covariates
NestHt <- as.vector(nest_covs$Nest_Ht)
InitDay <- as.vector(nest_covs$Init_day)
cosOrientation <- as.vector(nest_covs$cosOrientation) 
#Local covariates
Trees2550 <- as.vector(nest_covs$Trees_2550)
Trees50 <- as.vector(nest_covs$Trees_50)
PercPonderosa <- as.vector(nest_covs$pPIPO)
#climate covarites
PPT <- as.vector(nest_covs$meanPpt_mm)
Tmax <- as.vector(nest_covs$meanTmax_C)
#landscape-scale covariates
ForestCV <- as.vector(nest_covs$a1000_areacv2)
Contag <- as.vector(nest_covs$a1000_contag)
OpenNm <- as.vector(nest_covs$a1000_np1)
LandHa <- as.vector(nest_covs$a1000_Ha)
LandBu <- as.vector(nest_covs$a1000_RxBu)

# Exposure time -----------------------------------------------------------

#total amount of time each nest was surveyed
t <- as.vector(nests$exposure)

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
                 n.transects = n.transects,
                 n.trt = n.trt, 
                 n.species = n.species,
                 n.forests = n.forests,
                 n.stages = n.stages,
                 #Random effects IDs
                 Nest.num = Nest.num,
                 Year.num = Year.num,
                 Transect.num = Transect.num,
                 Forest.num = Forest.num,
                 #Missing data
                 Forest.ID = Forest.ID,
                 #Treatment covariate
                 TreatmentID = TreatmentID, 
                 #Nest species
                 SpeciesID = SpeciesID, 
                 #nest Stage
                 StageID = StageID,
                 #Nest-level covariates
                 NestHt = NestHt, 
                 cosOrientation = cosOrientation,
                 InitDay = InitDay, 
                 #Local-level covariates
                 Trees50 = Trees50,
                 Trees2550 = Trees2550, 
                 PercPonderosa = PercPonderosa,
                #climate covariates
                 Tmax = Tmax,
                 PPT = PPT,
                 #landscape-scale covariates
                 ForestCV = ForestCV,
                 Contag = Contag,
                 OpenNm = OpenNm,
                 LandHa = LandHa,
                 LandBu = LandBu,
                 #dataset
                 y = y,
                 #exposure
                 t = t)

saveRDS(all_data, here("data_outputs", 
                       "01_whwonests",
                       '03_JAGS_input_data',
                       "01_whwonests",
                      "mod1_JAGS_input_data.RDS"))

