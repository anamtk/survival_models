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
            "00_functions",
            "tidy_functions.R"))

# Load data ---------------------------------------------------------------

nests <- read.csv(here("data_outputs",
                       "01_whwonests", 
                       "02_analysis_ready",
                       "interval_models_nest_data.csv"))

full <- read.csv(here("data_outputs",
                       "01_whwonests",
                       "02_analysis_ready",
                       "full_survey_model_nest_data.csv"))

climate <- read.csv(here("data_outputs",
                         "01_whwonests",
                         "02_analysis_ready",
                         'daily_climate_data_withintervalID.csv'))


# Nest ID, number, and visit number ---------------------------------------

Nests <- nests %>%
  distinct(Nest_ID)


#get a vector of nest IDs
Nests <- unique(nests$Nest_ID)

# get the total count of nests
n.nests <- length(Nests)


#How many times did each nest get measured
#(number of intervals)
n.t <- nests %>%
  group_by(Nest_ID) %>%
  dplyr::select(Nest_ID, interval) %>%
  filter(interval == max(interval)) %>%
  #order these back in the right order bc they got confused
  arrange(ordered(Nest_ID, Nests)) %>%
  ungroup() %>%
  dplyr::select(interval) %>%
  as_vector()


mean(n.t)
sd(n.t)/sqrt(322)
max(n.t)
# Random variables of transect, forest, and year -------------------

#Year ID
Years <- nests %>%
  distinct(Year_located) %>%
  as_vector()

#number of years
n.years <- length(Years)

Forests <- nests %>% 
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

#year as random effect - vector length of nests
Year <- nests %>%
  distinct(Nest_ID, Year_located) %>%
  dplyr::select(Year_located) %>%
  as_vector()
#get numeric for model
Year.num <- nums(Year)

#forest as random effect - vector length of number 
# of transects because of hierarchical centering
Forest <- nests %>%
  distinct(Nest_ID, Project_ID) %>%
  dplyr::select(Project_ID) %>%
  as_vector()

#Numeric for model
Forest.num <- nums(Forest)

# Nest and stand covariates -----------------------------------------------
# **might be able to subset these based on previous literature**
#select all covariates on nest survival
nest_covs <- nests %>%
  distinct(Nest_ID, 
           Nest_Ht, Orientation, Initiation_date,
           Trees_2550, Trees_50) %>% 
  mutate(Initiation_date = case_when(str_detect(Initiation_date, "99") ~ NA_character_,
                                     TRUE ~ Initiation_date)) %>%
  mutate(Initiation_date = yday(Initiation_date)) %>%
  mutate_if(is.numeric, scale)  #center and scale continous variables

#Nest-level covariates
NestHt <- as.vector(nest_covs$Nest_Ht)
InitDay <- as.vector(nest_covs$Initiation_date)
cosOrientation <- as.vector(nest_covs$Orientation) 
#Local covariates
Trees2550 <- as.vector(nest_covs$Trees_2550)
Trees50 <- as.vector(nest_covs$Trees_50)


# Clmate data -------------------------------------------------------------


# Daily covariates --------------------------------------------------------

climate2 <- climate %>%
  arrange(match(Nest_ID, Nests)) %>%
  group_by(Nest_ID) %>%
  mutate(Survey_day = 1:n()) %>%
  ungroup()

n.days <- climate2 %>%
  group_by(Nest_ID) %>%
  tally() %>%
  arrange(match(Nest_ID, Nests)) %>%
  dplyr::select(n) %>%
  as_vector()


start.day <- climate2 %>%
  dplyr::select(Nest_ID, interval, Survey_day) %>%
  group_by(Nest_ID, interval) %>%
  filter(Survey_day == min(Survey_day)) %>%
  pivot_wider(names_from = 'interval',
              values_from = 'Survey_day') %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()

end.day <- climate2 %>%
  dplyr::select(Nest_ID, interval, Survey_day) %>%
  group_by(Nest_ID, interval) %>%
  filter(Survey_day == max(Survey_day)) %>%
  pivot_wider(names_from = 'interval',
              values_from = 'Survey_day') %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()

#[nest, day]
Tmax <- climate2 %>%
  dplyr::select(Nest_ID, tmax, Survey_day) %>%
  mutate(tmax = scale(tmax)) %>%
  pivot_wider(names_from = "Survey_day",
              values_from = "tmax")%>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()

PPT <- climate2 %>%
  dplyr::select(Nest_ID, ppt, Survey_day) %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = "Survey_day",
              values_from = "ppt")%>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()

# Calculate interval length for each interval
t <- nests %>%
  rowwise() %>%
  mutate(t = as_date(Visit_date) - as_date(Prev_date)) %>%
  mutate(t = as.numeric(t)) %>%
  dplyr::select(Nest_ID, interval, t) %>%
  pivot_wider(names_from = 'interval',
              values_from = "t") %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix()
  
# Response matrix ---------------------------------------------------------

#Observed interval success/failures for each nest
y <- nests %>%
  mutate(survival = case_when(Fate_cat == "success" ~ 1,
                              Fate_cat == "failure" ~ 0,
                              TRUE ~ NA_real_)) %>%
  group_by(Nest_ID) %>%
  dplyr::select(Nest_ID, survival, interval) %>%
  #fancy coding to set all previous visits to failed nests = 1
  #arrange so the last observation is first
  arrange(desc(interval)) %>%
  #get the last survival observation, and set all others = NA
  mutate(survival_last = case_when(duplicated(survival) ~ NA_real_,
                               TRUE ~ survival)) %>%
  #now re-create the survival column such that when it was the last observation
  # for that nest (survival_last not equal to NA), the nest gets the value for 
  # that last fate (so 1-0). If it was a nest observed in any but the last interval
  # the assumption is that the nest was alive - so those all get a value of 1
  mutate(survival = case_when(!is.na(survival_last) ~ survival_last, 
                               TRUE ~ 1)) %>%
  ungroup() %>%
  #arrange by interval
  arrange(interval) %>%
  #select variables for matrix
  dplyr::select(interval, survival, Nest_ID) %>%
  #make matrix
  pivot_wider(names_from = "interval",
              values_from = "survival") %>%
  #set nest names to rownames
  column_to_rownames(var = "Nest_ID") %>%
  #make a matrix
  as.matrix()
  
# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.nests = n.nests,
                 n.t = n.t, 
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
                 #daily data
                 n.days = n.days,
                 start.day = start.day,
                 end.day = end.day,
                 #interval lengths
                 t = t)

saveRDS(all_data, here("data_outputs", 
                       "01_whwonests",
                       '03_JAGS_input_data',
                       "mod2_daily_JAGS_input_data.RDS"))

