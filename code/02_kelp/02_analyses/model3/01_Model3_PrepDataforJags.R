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

kelp <- read.csv(here("data_outputs",
                      "02_kelp",
                      "02_analysis_ready",
                      "interval_kelp_data.csv"))

# Initial cleaning --------------------------------------------------------

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

# Arrange by interval number for custom model -----------------------

all_kelp2 <- all_kelp2 %>%
  group_by(Card_number) %>%
  mutate(n.t = n()) %>%
  ungroup() %>%
  arrange(n.t)


# Count variables for indexing --------------------------------------------

Plants <- all_kelp2 %>%
  distinct(Card_number, n.t)

#number with one interval only
n.plants1 <- all_kelp2 %>%
  filter(n.t == 1) %>%
  tally() %>%
  as_vector()

#total number of plants
n.plants <- length(Plants$Card_number)

#How many times did each nest get measured
#(number of intervals)
n.t <- as_vector(Plants$n.t)

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

n.substrates <- length(unique(all_kelp2$Substrate_Code))


# Random variables --------------------------------------------------------
#these next use the nums() function in the 
# custom data prep functions script sourced above

Transect.num <- all_kelp2 %>%
  unite(col = Site_transect, 
        c("Site", "Transect"),
        sep = "_" , 
        remove = FALSE) %>%
  distinct(Card_number, Site_transect) %>%
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


# Response data vector ----------------------------------------------------

y <- all_kelp2 %>%
  group_by(Card_number) %>%
  filter(Visit_interval == max(Visit_interval)) %>%
  ungroup() %>%
  distinct(Card_number, Presence) %>%
  dplyr::select(Presence) %>%
  as_vector()

length(y[which(y == 1)])
length(y[which(y ==0)])

# Export as RDS -----------------------------------------------------------

all_data <- list(n.plants1 = n.plants1, 
                 n.plants = n.plants,
                 n.t = n.t,
                 n.transects = n.transects,
                 n.sites = n.sites,
                 n.substrates = n.substrates,
                 Transect.num = Transect.num,
                 Site.num = Site.num,
                 SST = SST,
                 WavePower = WavePower,
                 Stipes = Stipes,
                 Diam = Diam,
                 Depth = Depth,
                 SubstrateID = SubstrateID,
                 t = t,
                 y = y)


saveRDS(all_data, here("data_outputs", 
                       "02_kelp",
                       '03_JAGS_input_data',
                       "mod3_JAGS_input_data.RDS"))

saveRDS(all_data, here("monsoon", 
                       "02_kelp",
                       "model3",
                       "inputs",
                       "mod3_JAGS_input_data.RDS"))

