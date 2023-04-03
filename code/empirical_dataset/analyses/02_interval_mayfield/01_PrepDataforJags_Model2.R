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

#nests8 dataset output from the cleaning code "04_nest_survival_dataprep.R"
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


nests1 <- nests1 %>%
  group_by(Nest_ID) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  arrange(total)
# Nest ID, number, and visit number ---------------------------------------

Nests <- nests1 %>%
  distinct(Nest_ID, total)

n.nests1 <- Nests %>%
  filter(total == 1) %>%
  tally() %>%
  as_vector()

n.nests2 <- Nests %>%
  filter(total > 1) %>%
  tally() %>%
  as_vector()

#get a vector of nest IDs
Nests <- unique(nests1$Nest_ID)

# get the total count of nests
n.nests <- length(Nests)


#How many times did each nest get measured
#(number of intervals)
n.t <- nests1 %>%
  group_by(Nest_ID) %>%
  tally(name = "n.t") %>%
  dplyr::select(n.t) %>%
  as_vector()

mean(n.t)
sd(n.t)/sqrt(318)
max(n.t)
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
           a1000_Bu, n_tx, Tm_since_tx, Time_groups) %>% 
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

TrtTime <- nest_covs %>%
  mutate(Time_groups = factor(Time_groups,
                             levels = c("oot", "0-3",
                                        "4-9", "10+"))) %>%
  dplyr::select(Time_groups) %>%
  as_vector() %>%
  nums()

n.times <- length(unique(TrtTime))

NTrt <- as.vector(nest_covs$n_tx)
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
LandBu <- as.vector(nest_covs$a1000_Bu)

# Sampling Interval Covariates --------------------------------------------

#all covariates collected at a per-interval period
JulianDate <- cov_matrix(df = nests1,  var = 'Julian_start')

#center and scale all the interval covariates
mJD <- mean(as.numeric(JulianDate), na.rm= T)
sdJD <- sd(as.numeric(JulianDate), na.rm = T)
JulianDate2 <- (JulianDate-mJD)/sdJD

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
  mutate(StageID = case_when(prevStage == "Ne" ~ 1,
                             prevStage == "Eg" ~ 2,
                             TRUE ~ NA_real_)) %>%
  group_by(Nest_ID) %>%
  mutate(interval  = row_number()) %>%
  ungroup() %>%
  dplyr::select(Nest_ID, StageID, interval) %>%
  pivot_wider(names_from = "interval",
              values_from = "StageID") %>%
  column_to_rownames(var = 'Nest_ID') %>%
  as.matrix()

n.stages <- length(unique(nests1$prevStage))
#get lag peeped value

nests1 %>%
  filter(Peeped == 1) %>%
  group_by(Nest_ID, Fate_cat) %>%
  tally() %>%
  ggplot(aes(x = Fate_cat, y = n)) +
  geom_boxplot()

#what this shows is that successful nests had more peeps,
# thus, misleading this estimate of survival based on
#whether the nest was peeped or not.

LagPeeped <- nests1 %>%
  mutate(survival = case_when(Fate_cat == "success" ~ 1,
                              Fate_cat == "failure" ~ 0,
                              TRUE ~ NA_real_)) %>%
  group_by(Nest_ID) %>%
  mutate(interval = row_number()) %>%
  dplyr::select(Nest_ID, Duration, Peeped, survival, interval,
                Year_located, Project_ID) %>%
  #fancy coding to set all previous visits to failed nests = 1
  #arrange so the last observation is first
  arrange(desc(interval)) %>%
  mutate(lag_peeped = lead(Peeped, n=1)) %>%
  arrange(interval) %>%
  ungroup() %>%
  dplyr::select(Nest_ID, interval, lag_peeped) %>%
  pivot_wider(names_from = "interval",
              values_from = "lag_peeped") %>%
  column_to_rownames(var = "Nest_ID") %>%
  mutate('1' = replace('1', is.na('1'), 0)) %>%
  mutate('1' = case_when('1' == 1 ~ 0,
                         TRUE ~ NA_real_)) %>%
  as.matrix() 

Age <- nests1 %>%
  rowwise() %>%
  mutate(Age = Julian_end - Init_day) %>%
  ungroup() %>%
  mutate(Age2 = scale(Age)) %>%
  group_by(Nest_ID) %>%
  mutate(interval = row_number()) %>%
  dplyr::select(Nest_ID, Age2,  interval) %>%
  pivot_wider(names_from = "interval",
              values_from = "Age2") %>%
  column_to_rownames(var = "Nest_ID") %>%
  as.matrix() 

# Response matrix ---------------------------------------------------------

#Observed interval success/failures for each nest
y <- nests1 %>%
  mutate(survival = case_when(Fate_cat == "success" ~ 1,
                              Fate_cat == "failure" ~ 0,
                              TRUE ~ NA_real_)) %>%
  group_by(Nest_ID) %>%
  mutate(interval = row_number()) %>%
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
  arrange(interval) %>%
  dplyr::select(interval, survival, Nest_ID) %>%
  pivot_wider(names_from = "interval",
              values_from = "survival") %>%
  column_to_rownames(var = "Nest_ID") 
  
# Variables for derived survival estimates --------------------------------

#how many days are nests in each group? 
# these values are taken from the literature
LayDays <- 5
IncubDays <- 14
YoungDays <- 26
#around 45 day total nest period

#Median Nesting day for prob calculations
nests1 %>%
  summarise(median = median(Julian_start, na.rm = T))
#177
# full period is 40 days according to the literature
# But 35 according to our data:
nests1 %>%
  group_by(Stage) %>%
  summarise(mean = mean(Julian_start, na.rm = T))

# 177 is median of the nesting period, so
# nesting period ranges from 155 - 192
#To get the median Julian day for each period, 
#get the median of the day range for each period
# and add the start date
median(c(1,2,3)) + 155#2
median(c(4:19)) +155 #11.5
median(c(20:40)) + 155 #30

#then we can center those
Lay.Med.JD <- (157 - mJD)/sdJD
Incub.Med.JD <- (166.5-mJD)/sdJD
Young.Med.JD <- (185-mJD)/sdJD

# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.nests = n.nests,
                 n.t = n.t, 
                 n.years = n.years,
                 n.transects = n.transects,
                 n.trt = n.trt, 
                 n.species = n.species,
                 n.stages = n.stages,
                 n.forests = n.forests,
                 n.times = n.times,
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
                 TrtTime = TrtTime,
                 NTrt = NTrt,
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
                 #interval lengths
                 t = t)

saveRDS(all_data, here("data_outputs", 
                       'model_input_data',
                      "survival_JAGS_input_data.RDS"))


# Data summaries ----------------------------------------------------------
#how many nests per year:
nests1 %>%
  distinct(Nest_ID, Year_located, Project_ID, Trt_cat, Fate,
           Fate_cat) %>%
  group_by(Year_located) %>%
  tally()

#avg per year
nests1 %>%
  distinct(Nest_ID, Year_located, Project_ID, Trt_cat, Fate,
           Fate_cat) %>%
  group_by(Year_located) %>%
  tally() %>%
  summarise(mean = mean(n),
            sd = sd(n),
            total = n(),
            se = sd/sqrt(total))

#how many nests per treatment category? 
nests1 %>%
  distinct(Nest_ID, Year_located, Project_ID, Trt_cat, Fate,
           Fate_cat) %>%
  group_by(Trt_cat) %>%
  tally() 

#survival overall
nests1 %>%
  distinct(Nest_ID, Year_located, Project_ID, Trt_cat, Fate,
           Fate_cat) %>%
  group_by(Fate_cat) %>%
  tally()

#survival by year 
nests1 %>%
  distinct(Nest_ID, Year_located, Project_ID, Trt_cat, Fate,
           Fate_cat) %>%
  group_by(Fate_cat, Year_located) %>%
  tally() %>%
  ungroup() %>%
  pivot_wider(names_from = "Fate_cat",
              values_from = "n") %>%
  rowwise() %>%
  summarise(success_rate = success/(success + failure)) %>%
  summarise(mean_rate = mean(success_rate,na.rm =T),
            sd = sd(success_rate, na.rm =T),
            total = n(),
            se = sd/sqrt(total))

# cause of mortality
nests1 %>%
  distinct(Nest_ID, Year_located, Project_ID, Trt_cat, Fate,
           Fate_cat) %>%
  group_by(Fate) %>%
  tally()
#1 = success
#2 = bear
# 3 = corvid
# 4 = squirrel
# 5 = chipmunk
# 6 = snake
# 7 = weather
# 8 = cavity destroyed
# 9 = unknown
# 10 = other

