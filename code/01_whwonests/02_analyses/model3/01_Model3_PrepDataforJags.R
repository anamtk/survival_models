# Nest survival data prep
# November 1, 2021
# Ana Miller-ter Kuile

# prepping data for the model of nest survival

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "lubridate")

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

climate <- read.csv(here("data_outputs",
                         "01_whwonests",
                         "02_analysis_ready",
                         'daily_climate_data_withintervalID.csv'))

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

#total number of nests
n.nests <- length(Nests$Nest_ID)

#How many times did each nest get measured
#(number of intervals)
n.t <- as_vector(Nests$n.t)

# Random variables of transect, forest, and year -------------------

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

#year as random effect - vector length of nests
Year <- nests1 %>%
  distinct(Nest_ID, Year_located) %>%
  dplyr::select(Year_located) %>%
  as_vector()
#make numeric for model
Year.num <- nums(Year)

#forest ID length of number of transects for hierarchical centered
Forest <- nests1 %>%
  distinct(Nest_ID,
           Project_ID) %>%
  dplyr::select(Project_ID) %>%
  as_vector()
#make numeric for model
Forest.num <- nums(Forest)

# Nest and stand covariates -----------------------------------------------
# **might be able to subset these based on previous literature**
#select all covariates on nest survival
nest_covs <- nests1 %>%
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


# Daily covariates --------------------------------------------------------

climate2 <- climate %>%
  arrange(match(Nest_ID, Nests$Nest_ID)) %>%
  group_by(Nest_ID) %>%
  mutate(Survey_day = 1:n()) %>%
  ungroup()

n.days <- climate2 %>%
  group_by(Nest_ID) %>%
  tally() %>%
  arrange(match(Nest_ID, Nests$Nest_ID)) %>%
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

# Response vector ---------------------------------------------------------

y <- nests1 %>%
  distinct(Nest_ID, Fate_cat) %>%
  mutate(survival = case_when(Fate_cat == "success" ~ 1,
                              Fate_cat == "failure" ~ 0,
                              TRUE ~ NA_real_)) #%>%
  dplyr::select(survival) %>% 
  as_vector()
  
length(y[which(y == 1)]) #238
length(y[which(y == 0)]) #84

# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.nests1 = n.nests1,
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
                 #duplicated dataset for ifelse statement
                 y2 = y,
                 #interval and day indexing
                 start.day = start.day,
                 end.day = end.day,
                 n.days = n.days)

saveRDS(all_data, here("data_outputs", 
                       "01_whwonests",
                       '03_JAGS_input_data',
                       "mod3_JAGS_input_data.RDS"))

