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

kelp <- read.csv(here("data_outputs",
                      "02_kelp",
                      "02_analysis_ready",
                       "interval_kelp_data.csv"))

# Variables of interest ---------------------------------------------------

#a couple plants only visited once
all_kelp2 <- kelp %>%
  group_by(Card_number) %>%
  mutate(sum = n()) %>%
  filter(sum > 1) %>%
  ungroup() %>%
  dplyr::select(-sum) %>%
  #also - remove first visit interval since all alive then
  filter(Visit_interval != 0)  #3799

#1838
all_kelp2 %>%
  distinct(Card_number) %>%
  tally()

#2 plants removed that were only visited once
kelp %>%
  distinct(Card_number) %>%
  tally()

# Count variables for indexing --------------------------------------------


#count variables for indexing
n.plants <- all_kelp2 %>%
  distinct(Card_number) %>%
  summarise(n.plants = n()) %>%
  as_vector()

n.transects <- all_kelp2 %>%
  unite(col = Site_transect, 
        c("Site", "Transect"),
        sep = "_" , 
        remove = FALSE) %>%
  distinct(Site_transect) %>%
  summarise(n.transects = n()) %>%
  as_vector()

n.sites <- all_kelp2 %>%
  distinct(Site) %>%
  summarise(n.reefs = n()) %>%
  as_vector()

#this should be a by-plant vector of number of visit intervals
n.t <- all_kelp2 %>%
  mutate(Card_number = factor(Card_number, levels = unique(Card_number))) %>%
  group_by(Card_number) %>%
  summarise(n.t = n()) %>%
  dplyr::select(n.t) %>%
  as_vector()

n.substrates <- length(unique(all_kelp2$Substrate_Code))

# Random variables --------------------------------------------------------
#these next use the nums() function in the 
# custom data prep functions script sourced above

Transect.num <- all_kelp2 %>%
  unite(col = Site_transect, 
        c("Site", "Transect"),
        sep = "_" , 
        remove = FALSE) %>%
  dplyr::select(Site_transect) %>%
  as_vector() %>%
  nums()

#site num length of number of transects - 15
Site.num <- all_kelp2 %>%
  unite(col = Site_transect, 
        c("Site", "Transect"),
        sep = "_" , 
        remove = FALSE) %>%
  distinct(Site, Site_transect) %>%
  dplyr::select(Site) %>%
  as_vector() %>%
  nums()

# Interval covariates -----------------------------------------------------
#get a scaled matrix for each interval for each plant
#these use the int_cov custom function from the 
#data prep functions R script

SST <- int_cov(variable = sst_mean)

WavePower <- int_cov(variable = wave_p_mean)

Stipes <- int_cov(variable = prev_Stipes)

# Fixed covariates for each plant -----------------------------------------

#these are fixed covariates for each plant
# and derived via the fixed_covs function
# fom the data prep functions R script
Diam <- fixed_covs(variable = Max_Diam_1_cm)

Depth <- fixed_covs(variable = Depth_ft)

#Substrate
SubstrateID <- all_kelp2 %>%
  distinct(Card_number, Substrate_Code) %>%
  mutate(Substrate_Code = as.numeric(as.factor(Substrate_Code))) %>%
  dplyr::select(Substrate_Code) %>%
  as_vector()

#check that substrate code with most observations is
# the baseline code
all_kelp2 %>%
  distinct(Card_number, Substrate_Code) %>%
  group_by(Substrate_Code) %>%
  tally()

levels(as.factor(all_kelp2$Substrate_Code))
#B = 1 = bedrock
#BO = 2 = boulders
#C = 3 = cobble
#S = 4 = sand
#SS = 5 = shallow sand

hist(SubstrateID) #lots of 1, 6, 7 - may be a problem 
all_kelp2 %>%
  distinct(Card_number, Substrate_Code) %>%
  group_by(Substrate_Code) %>%
  tally()

# Time interval matrix ----------------------------------------------------

# this should be a plant x interval matrix
t <- all_kelp2 %>%
  dplyr::select(Card_number, 
                Visit_Date,
                prev_Visit_Date, 
                Visit_interval) %>%
  mutate(interval_length = as.Date(Visit_Date) - 
           as.Date(prev_Visit_Date)) %>%
  mutate(interval_length = as.numeric(interval_length)) %>%
  dplyr::select(Card_number, Visit_interval, 
                interval_length) %>%
  pivot_wider(names_from = Visit_interval,
              values_from = interval_length) %>%
  column_to_rownames(var = 'Card_number') %>%
  as.matrix()


# Response data y matrix --------------------------------------------------

#underlyin data should be a plant x interval matrix
y <- all_kelp2 %>%
  dplyr::select(Card_number, Visit_interval,
                Presence) %>%
  pivot_wider(names_from = Visit_interval,
              values_from = Presence) %>%
  column_to_rownames(var = "Card_number") %>%
  as.matrix()


# Compile and export ------------------------------------------------------

all_data <- list(n.plants = n.plants,
                 n.t = n.t,
                 n.transects = n.transects,
                 n.substrates = n.substrates,
                 n.sites = n.sites,
                 Transect.num = Transect.num,
                 Site.num = Site.num,
                 SST = SST, 
                 Stipes = Stipes,
                 WavePower = WavePower,
                 Diam = Diam,
                 Depth = Depth,
                 SubstrateID = SubstrateID,
                 t = t,
                 y = y)


saveRDS(all_data, here("data_outputs", 
                       "02_kelp",
                       "03_JAGS_input_data",
                       "mod2_JAGS_input_data.RDS"))

saveRDS(all_data, 
        file = here("monsoon",
                    "02_kelp",
                    "model2",
                    "inputs",
                    "mod2_JAGS_input_data.RDS"))





