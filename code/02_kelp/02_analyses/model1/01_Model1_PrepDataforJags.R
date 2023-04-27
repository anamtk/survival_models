# Nest survival data prep - full interval model1
# Ana Miller-ter Kuile
# March 31, 2023

# prepping data for the model of nest survival where data are
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

kelp <- read.csv(here("data_outputs",
                      "02_kelp",
                      "02_analysis_ready",
                      "total_intervals_kelp_data.csv"))

# Plant ID, number --------------------------------------------------------

Kelp <- kelp %>%
  distinct(Card_number)

#get a vector of plant IDS
Kelp <- unique(kelp$Card_number)

# get the total count of plants
n.plants <- length(Kelp)

# Random variables of transect, forest ------------------------------------

#Get the names of all the transects
Transects <- kelp %>%
  distinct(Site_transect) %>%
  as_vector()

#number of transects
n.transects <- length(Transects)

Sites <- kelp %>%
  distinct(Site) %>%
  as_vector()

#number of sites
n.sites <- length(Sites)

# Create covariates -------------------------------------------------------

#the next sections create variables that encompass the random
# variables that group the data,
# and covariates that relate to the environment

# Random variables --------------------------------------------------------

#transect random effect (vector length of number of plants)
Transect <- kelp %>%
  distinct(Card_number,
           Site_transect,
           Site) %>%
  dplyr::select(Site_transect) %>%
  as_vector()

#make numeric for the model
Transect.num <- nums(Transect)

#site random effect - vector of length of the 
#number of transects (hierarchically centered)
Site <- kelp %>%
  distinct(Site_transect,
           Site) %>%
  dplyr::select(Site) %>%
  as_vector()

#make numeric for model
Site.num <- nums(Site)


# Covariates --------------------------------------------------------------

covs <- kelp %>%
  mutate(Card_number = as.factor(Card_number)) %>%
  distinct(Card_number, Depth_ft,
                Max_Diam_1_cm, Stipe_num,
                sst_mean, wave_p_mean,
                Substrate_Code) %>%
  mutate(Depth_m = Depth_ft*0.3048) %>%
  dplyr::select(-Depth_ft) %>%
  mutate_if(is.numeric, scale) #center and scale continuous covariates
  

SST <- as.vector(covs$sst_mean)
WavePower <- as.vector(covs$wave_p_mean)
Stipes <- as.vector(covs$Stipe_num)
Diam <- as.vector(covs$Max_Diam_1_cm)
Depth <- as.vector(covs$Depth_m)

covs %>%
  group_by(Substrate_Code) %>%
  tally()

SubstrateID <- covs %>%
  mutate(Substrate_Code = factor(Substrate_Code,
                                 levels = c("B", "BO",
                                            "C", "S", "SS"
                                            ))) %>%
  dplyr::select(Substrate_Code) %>%
  as_vector() %>%
  nums()

n.substrates <- length(unique(SubstrateID))

# Exposure time -----------------------------------------------------------

#total amount of time each nest was surveyed
t <- as.vector(kelp$t)

# Response vector ---------------------------------------------------------

y <- kelp  %>%
  dplyr::select(Presence) %>%
  #make a vector
  as_vector()

# Export as RDS -----------------------------------------------------------

all_data <- list(#Data count variables
                 n.plants = n.plants,
                 n.transects = n.transects,
                 n.sites = n.sites,
                 n.substrates = n.substrates,
                 #Random effects IDs
                 Transect.num = Transect.num,
                 Site.num = Site.num,
                 #Covariates
                 SST = SST,
                 WavePower = WavePower,
                 Stipes = Stipes,
                 Diam = Diam,
                 Depth = Depth,
                 SubstrateID = SubstrateID,
                 #dataset
                 y = y,
                 #exposure
                 t = t)

saveRDS(all_data, here("data_outputs", 
                       "02_kelp",
                       "03_JAGS_input_data",
                      "mod1_JAGS_input_data.RDS"))

