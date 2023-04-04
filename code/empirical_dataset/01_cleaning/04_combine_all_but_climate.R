# Nest survival data prep
# Ana Miller-ter Kuile
# October 28, 2021

# this script prepares data from nests where fate 
# was tracked to do survival analysis, exploring
# factors that influence the survival of nests.

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "readxl", "lubridate")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#temp and precip functions
source(here("code",
            "functions",
            "tidy_functions.R"))

# Load data ---------------------------------------------------------------

#get the dataframe with all the treatment info for each nest
trt_covs <- read.csv(here("data_outputs",
                          "01_cleaning",
                          "01_Nest_trt_covs.csv"))

#tree data compilation, including abundance and sizes of trees and snags
tree_data <- read.csv(here("data_outputs",
                                "01_cleaning",
                                "02_nest_tree_densities.csv"))


#landscape variables
land <- read.csv(here("data_outputs",
                      "01_cleaning",
                      "03_nest_landscape_1000_covariates.csv"))

#nest locations DF
fates <- read_xlsx(here("data_raw",
                        "bird_data",
                        "Birds03_nest_fates.xlsx"))

#Individual nest visits DF
visits <- read_xlsx(here("data_raw",
                         "bird_data",
                         "Birds02_nest_visits.xlsx"))
 

# Combine nest location,tx, and tree data -------------------------------------

nests <- trt_covs %>%
  dplyr::select(-Point_ID) %>%
  left_join(tree_data, by = c("Nest_ID" = "Measurement_ID")) 

#how many nests are in this df?
nests %>%
  tally()
#363

# Combine fates with nest and treatment data ------------------------------

nests2 <- nests %>%
  left_join(fates, by = c("Nest_ID", 'Date_located')) %>%
  #cateogrize into success/failure, and unknown
  mutate(Fate_cat = case_when(Fate == 1 ~ "success",
                              Fate == 11 ~ "unknown",
                              TRUE ~ "failure")) %>%
  #remove unknown for future analyses
  filter(Fate_cat != "unknown") 
#This may remove nests that had initial egg counts in them -
#I wonder if it is worth keeping those in?? not sure

nests2 %>%
  tally()
#330 nests
# Clean and combine with visits dataset -----------------------------------

#get nest IDs so we can filter the visits dataset for just nests
nest_IDs <- nests2$Nest_ID

#Get visits dataset to be just nest visits
visits2 <- visits %>%
  # find visits that correspond to a Nest ID
  filter(Nest_ID %in% c(nest_IDs)) %>%
  #select variables of interest
  dplyr::select(Nest_ID, Visit_date, No_eggs, No_young, Stage, 
                Peeped, St_time, End_time) %>%
  #figure out the time spent at each nest during each visit
  mutate(interval = St_time %--% End_time,
         #convert to a duration in minutes
         Duration = as.duration(interval)/dminutes(1)) %>%
  # fix typos to be NA values for now 
  mutate(Duration = case_when(Duration < 0 ~ NA_real_, 
                              TRUE ~ Duration))

# Select variables of interest, convert to Julian days and tiem si --------

nests3 <- visits2 %>%
  #join the visits dataset with all the nest metatdata from nests
  left_join(nests2, by = "Nest_ID") %>%
  #select variables of interest
  dplyr::select(Nest_ID, Visit_date, No_eggs, No_young,
                Stage, Year_located, Project_ID, 
                Transect_ID, Trt_cat, Fate,
                UTM_datum_zone, UTM_N, UTM_E, Nest_Ht,
                Decay_class, Tree_Snag_Log,
                Tree_sp, Tree_ht, DBH, Orientation, 
                Aspect,
                Slope, Initiation_date, 
                NoFL_uncert, NoFL_cert, Trees_2550,  Trees_50,
                pPIPO, PIPO_2550, PIPO_50,  Fate_cat,
                Peeped, Duration) %>%
  # get Julian Days for each start and end
  mutate(Julian_end = as.Date(Visit_date)) %>%
  mutate(Julian_end = format(Julian_end, "%j")) %>%
  group_by(Nest_ID) %>%
  #arrange visits in order by Julian Start to calculate end (from next visit)
  arrange(Julian_end, .by_group = TRUE) %>%
  #get end date based on the next visit period
  mutate(Julian_start = lag(Julian_end, 1)) %>%
  mutate(Initiation_date = ymd(Initiation_date)) %>%
  mutate(Init_day = yday(Initiation_date)) %>%
  ungroup() %>%
  #remove unknown stage nests 
  filter(!Stage %in% c("9")) %>%
  mutate(Orientation = case_when(Orientation == 999 ~ NA_integer_, 
                                 TRUE ~ Orientation)) %>%
  #set orientation to be 1 if N, -1 if South
  mutate(cosOrientation = cos(Orientation * (pi/180)))

# t <- nests3 %>%
#   dplyr::select(Nest_ID, Julian_start, Julian_end, Init_day,
#                 Visit_date, Initiation_date, Stage)
# 
# t2 <- t %>%
#   filter(Julian_end < Init_day)
# Remove multiple visits per day ------------------------------------------

#some nests were visited twice in one day and the number of eggs/young
# counted on the second visit. We want ot delete the duplicate
# with the missing info

nests4 <- nests3 %>%
  group_by(Nest_ID, Visit_date) %>%
  #get nests where visits are greater than 1
  mutate(sum = n()) %>%
  # then delete the ones where they were visited twice ANd eggs weren't counted
  filter(!(sum > 1 & No_eggs == 999)) %>%
  #filter(Trt_150 != "B") #remove one nest where treatment is burn
  ungroup()

# final total of 330 nests
nests4 %>%
  distinct(Nest_ID) %>%
  tally()

# Combine landscape variables ---------------------------------------------

nests5 <- nests4 %>%
  left_join(land, by = "Nest_ID")

# Recategorize nest tree IDs ----------------------------------------------
#POTR5 - Populus tremuloides
#PIPO - Ponderosa pine
#JUOC - Juniper occidentalis
#PSME/PSMEG = doug fir
#ABGR/ABCO = firs, depending on location

#PSMEG and PSME ==
#ABCO and AB... ==

nests5 <- nests5 %>%
  mutate(Tree_sp = case_when(Tree_sp %in% c("PSMEG", "PSME") ~ "PSME",
                             Tree_sp %in% c("ABCO", "ABGR") ~ "Abies",
                             TRUE ~ Tree_sp))



# Add a new transect ID for incidental transects --------------------------

#in order for the model to not lump all "incidental" 
# nest spatially with the hierarchical effect - setting
#each "incidental" nest it's own "transectID"

transect_ids <- nests5 %>%
  distinct(Nest_ID, Transect_ID) %>%
  mutate(num = as.character(1:n())) %>%
  mutate(Transect_ID2 = case_when(Transect_ID == "Incidental" ~ num,
                                  TRUE ~ Transect_ID)) %>%
  dplyr:: select(-num, - Transect_ID)

nests5 <- nests5 %>%
  left_join(transect_ids, by = "Nest_ID")


# Clean up columns in datafram --------------------------------------------

colnames(nests5)

nests6 <- nests5 %>%
  dplyr::select(Nest_ID, Visit_date,No_eggs, No_young,
                Stage, Year_located, Project_ID, Transect_ID,
                Transect_ID2,
                Trt_cat, Fate, UTM_datum_zone, UTM_N, UTM_E,
                Nest_Ht, Tree_sp, Orientation, cosOrientation,
                Initiation_date,Init_day,
                NoFL_uncert, NoFL_cert, Trees_2550, Trees_50,
                pPIPO, Fate_cat, Peeped, Duration, Julian_end,
                Julian_start, 
                a1000_areamn2, a1000_areaam2,
                a1000_areacv2, a1000_contag,
                a1000_pland2,
                a1000_pland1,
                a1000_lpi, a1000_lpi1, a1000_lpi2,
                a1000_np, a1000_np1, a1000_np2,
                a1000_proxmn2,a1000_Ha, a1000_RxBu,
                a1000_WFBu, a1000_PBu)

#filter out wildfire nests
nests6 <- nests6 %>%
  filter(Trt_cat != "W")

# Double check all variables of interest ----------------------------------

data_check <- nests6 %>%
  dplyr::select(Nest_Ht,Orientation, Init_day,
                Trees_2550, Trees_50, pPIPO,
                a1000_areamn2:a1000_RxBu
                )


ggplot(gather(data_check), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free')
# Export ------------------------------------------------------------------

write.csv(nests6, 
          here('data_outputs', 
               "01_cleaning", 
               "04_nest_local_land_data.csv"))


# Export for Kyle climate -------------------------------------------------

colnames(nests6)

data <- nests6 %>%
  dplyr::select(Nest_ID, 
                Year_located, Julian_start, Julian_end, 
                UTM_datum_zone, UTM_N, UTM_E)

write.csv(data, 
          here('data_outputs', 
               "01_cleaning", 
               "04a_nest_locations_dates.csv"))

# END SCRIPT ##

