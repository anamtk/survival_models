# Nest survival data prep
# November 1, 2021
# Ana Miller-ter Kuile

# prepping data for the model of nest survival

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

nests <- read.csv(here("data_outputs",
                       "02_analysis_ready",
                       "empirical", 
                       "interval_models_nest_data.csv"))


# Arrange nests by interval number for custom model -----------------------

nests1 <- nests %>%
  group_by(Nest_ID) %>%
  mutate(n.t = n()) %>%
  ungroup() %>%
  arrange(n.t)
# Nest ID, number, and visit number ---------------------------------------

Nests <- nests1 %>%
  distinct(Nest_ID, n.t)

#number with one interval only
n.nests1 <- Nests %>%
  filter(n.t == 1) %>%
  tally() %>%
  as_vector()

#number with more than one interval
n.nests2 <- Nests %>%
  filter(n.t > 1) %>%
  tally() %>%
  as_vector()

tot.nests <- length(Nests$Nest_ID)

#How many times did each nest get measured
#(number of intervals)
n.t <- as_vector(nests1$n.t)

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

Nest.num <- 1:tot.nests

Forest <- nests1 %>%
  distinct(Nest_ID,
           Transect_ID2,
           Project_ID) %>%
  dplyr::select(Project_ID) %>%
  as_vector()
Forest.num <- nums(Forest)
           
# Nest and stand covariates -----------------------------------------------
# **might be able to subset these based on previous literature**
#select all covariates on nest survival
nest_covs <- nests1 %>%
  distinct(Nest_ID, 
           Nest_Ht, Tree_sp, cosOrientation, Init_day,
          Trt_cat,
           Trees_2550, Trees_50, pPIPO, 
           a1000_areacv2, a1000_contag,
           a1000_np1, a1000_Ha, a1000_RxBu) %>% 
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
ForestCV <- as.vector(nest_covs$a1000_areacv2)
Contag <- as.vector(nest_covs$a1000_contag)
OpenNm <- as.vector(nest_covs$a1000_np1)
LandHa <- as.vector(nest_covs$a1000_Ha)
LandBu <- as.vector(nest_covs$a1000_RxBu)

# Sampling Interval Covariates --------------------------------------------

#Temperature and precip
Tmax <- nests1 %>%
  mutate(meanTmax_C = scale(meanTmax_C)) %>%
  group_by(Nest_ID) %>%
  mutate(interval = row_number()) %>%
  ungroup() %>%
  dplyr::select(Nest_ID, meanTmax_C, interval) %>%
  pivot_wider(names_from= "interval",
              values_from = "meanTmax_C") %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()


PPT <- nests1 %>%
  mutate(meanPpt_mm = scale(meanPpt_mm)) %>%
  group_by(Nest_ID) %>%
  mutate(interval = row_number()) %>%
  ungroup() %>%
  dplyr::select(Nest_ID, meanPpt_mm, interval) %>%
  pivot_wider(names_from= "interval",
              values_from = "meanPpt_mm") %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()

# Calculate interval length for each interval
t <- nests1 %>%
  mutate(Julian_end = as.numeric(Julian_end),
         Julian_start = as.numeric(Julian_start)) %>%
  mutate(Int = Julian_end - Julian_start) %>%
  group_by(Nest_ID) %>%
  mutate(interval = row_number()) %>%
  ungroup() %>%
  dplyr::select(Nest_ID, interval, Int) %>%
  pivot_wider(names_from = "interval",
              values_from = "Int") %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()
  

nests1 %>%
  group_by(prevStage) %>%
  tally()

#1 = Ne = 1186
#2 = Eg = 381
#3 = Ex = 108
# get stage of nest for each visit too
#this should be previous stage
Stage <- nests1 %>%
  mutate(StageID = case_when(prevStage == "N" ~ 1,
                             prevStage %in% c("I", "L") ~ 2,
                             TRUE ~ NA_real_)) %>%
  group_by(Nest_ID) %>%
  mutate(interval  = row_number()) %>%
  ungroup() %>%
  dplyr::select(Nest_ID, StageID, interval) %>%
  pivot_wider(names_from = "interval",
              values_from = "StageID") %>%
  column_to_rownames(var = 'Nest_ID') %>%
  as.matrix()

n.stages <- 2

Age <- nests1 %>%
  mutate(Age2 = scale(Age)) %>%
  group_by(Nest_ID) %>%
  mutate(interval = row_number()) %>%
  dplyr::select(Nest_ID, Age2,  interval) %>%
  pivot_wider(names_from = "interval",
              values_from = "Age2") %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix() 


# Response vector ---------------------------------------------------------

y <- nests1 %>%
  distinct(Nest_ID, Fate_cat) %>%
  mutate(survival = case_when(Fate_cat == "success" ~ 1,
                              Fate_cat == "failure" ~ 0,
                              TRUE ~ NA_real_)) %>%
  dplyr::select(survival) %>% 
  as_vector()
  
length(y[which(y == 1)]) #228
length(y[which(y == 0)]) #92

# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.nests1 = n.nests1,
                 n.nests2 = n.nests2,
                 tot.nests = tot.nests,
                 n.t = n.t, 
                 n.years = n.years,
                 n.transects = n.transects,
                 n.trt = n.trt, 
                 n.species = n.species,
                 n.stages = n.stages,
                 n.forests = n.forests,
                 #Random effects IDs
                 Nest.num = Nest.num,
                 Year.num = Year.num,
                 Transect.num = Transect.num,
                 Forest.num = Forest.num,
                 #Interval covariate
                 StageID = Stage,
                 Age = Age,
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
                 y = y,
                 y2 = y,
                 #interval lengths
                 t = t)

saveRDS(all_data, here("data_outputs", 
                       '03_JAGS_input_data',
                       "empirical",
                       "mod3_JAGS_input_data.RDS"))

