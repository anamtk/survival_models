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
            "functions",
            "tidy_functions.R"))

# Load data ---------------------------------------------------------------

nests <- read.csv(here("data", 
                       "Nest_survival_data.csv"))


# Remove specific observations --------------------------------------------

nests1 <- nests %>%
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

#318 nests total now
nests1 %>%
  distinct(Nest_ID) %>%
  tally()

# Nest ID, number, and visit number ---------------------------------------

Nests <- nests1 %>%
  distinct(Nest_ID, total)

#get a vector of nest IDs
Nests <- unique(nests1$Nest_ID)

# get the total count of nests
n.nests <- length(Nests)

# Random variables of transect, forest, and year -------------------

#Get the names of all the transects
Transects <- nests1 %>%
  distinct(Transect_ID2) %>%
  as_vector()

#number of transects
n.transects <- length(Transects)

#Year ID
Years <- nests1 %>%
  distinct(Year_located) %>%
  as_vector()

#number of years
n.years <- length(Years)

Forests <- nests1 %>% 
  distinct(Project_ID) %>%
  as_vector()

n.forests <- length(Forests)

# Create covariates -------------------------------------------------------

#the next sections create covariates that encompass the random
# variables that group the data,
# covariates that are at the different spatial scales 
# (e.g. tree density, nest location variables, forest treatment variables)
# and covariates that are at the sampling interval level
# (e.g. how long was the interval,
# what stage was the nest in during that interval)

# Random variables --------------------------------------------------------

#transect random effect (vector length of number of nest points)
Transect <- nests1 %>%
  distinct(Nest_ID, 
           Transect_ID2,
           Project_ID) %>%
  dplyr::select(Transect_ID2) %>%
  as_vector()
Transect.num <- nums(Transect)


#year as random effect - vector length of nests
Year <- nests1 %>%
  distinct(Nest_ID, Year_located) %>%
  dplyr::select(Year_located) %>%
  as_vector()
Year.num <- nums(Year)

Nest.num <- 1:n.nests

Forest <- nests1 %>%
  distinct(Nest_ID,
           Transect_ID,
           Project_ID) %>%
  dplyr::select(Project_ID) %>%
  as_vector()
Forest.num <- nums(Forest)
           
# Nest and stand covariates -----------------------------------------------
# **might be able to subset these based on previous literature**
#select all covariates on nest survival
nest_covs <- nests1 %>%
  distinct(Nest_ID, Nest_Ht, 
           Tree_sp, Orientation, Init_day,
          pPIPO, Trt_cat,
           Trees_2550, Trees_50,
           Project_ID, 
           Tmax, PPT,
           a1000_areacv2, a1000_contag,
           a1000_np1, a1000_Ha,
           a1000_RxBu, n_tx, Tm_since_tx, Time_groups) %>% 
  mutate(cosOrientation = cos(Orientation*(pi/180))) %>% #1 = north, -1 = south
  mutate_if(is.numeric, scale)  #center and scale continous variables

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

#Nest-level covariates
NestHt <- as.vector(nest_covs$Nest_Ht)

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

InitDay <- as.vector(nest_covs$Init_day)
cosOrientation <- as.vector(nest_covs$cosOrientation) 

#Local covariates
Trees2550 <- as.vector(nest_covs$Trees_2550)
Trees50 <- as.vector(nest_covs$Trees_50)
PercPonderosa <- as.vector(nest_covs$pPIPO)

#landscape-scale covariates
PPT <- as.vector(nest_covs$PPT)
Tmax <- as.vector(nest_covs$Tmax)
ForestCV <- as.vector(nest_covs$a1000_areacv2)
Contag <- as.vector(nest_covs$a1000_contag)
OpenNm <- as.vector(nest_covs$a1000_np1)
LandHa <- as.vector(nest_covs$a1000_Ha)
LandBu <- as.vector(nest_covs$a1000_RxBu)


# Response vector ---------------------------------------------------------

y <- nests1  %>%
  distinct(Nest_ID, Fate_cat) %>%
  mutate(survival = case_when(Fate_cat == "success" ~ 1,
                              Fate_cat == "failure" ~ 0,
                              TRUE ~ NA_real_)) %>%
  dplyr::select(survival) %>%
  as_vector()

# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.nests = n.nests,
                 n.years = n.years,
                 n.transects = n.transects,
                 n.trt = n.trt, 
                 n.species = n.species,
                 n.forests = n.forests,
                 #Random effects IDs
                 Nest.num = Nest.num,
                 Year.num = Year.num,
                 Transect.num = Transect.num,
                 Forest.num = Forest.num,
                 #Treatment covariate
                 TreatmentID = TreatmentID, 
                 #Nest-level covariates
                 NestHt = NestHt, 
                 cosOrientation = cosOrientation,
                 InitDay = InitDay, 
                 SpeciesID = SpeciesID, 
                 #Local-level covariates
                 Trees50 = Trees50,
                 Trees2550 = Trees2550, 
                 PercPonderosa = PercPonderosa,
                 #landscape-scale covariates
                 Tmax = Tmax,
                 PPT = PPT,
                 ForestCV = ForestCV,
                 Contag = Contag,
                 OpenNm = OpenNm,
                 LandHa = LandHa,
                 LandBu = LandBu,
                 #dataset
                 y = y)

saveRDS(all_data, here("data_outputs", 
                       'model_input_data',
                       "empirical",
                       "model1",
                      "survival_mod1_JAGS_input_data.RDS"))

