# Nest data compiling
# March 2, 2022
# Ana Miller-ter Kuile

# this script uses QCed nest data with associated treatment codes
# to generate the treatment and nest-level covariates for the 
# nest model


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse", 
                  "readxl")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

#need: locations dataframe

locations <- read_xlsx(here("data_raw",
                            "bird_data",
                            "Birds01_nest_locations.xlsx"))


# FACTS dataset
facts <- read_xlsx(here("data_raw",
                        "bird_data",
                        "Birds07_nest_FACTS.xlsx"))


# Select CFLRP nests ------------------------------------------------------

#select locations only in the three CFLRP programs
locs <- locations %>%
  filter(Project_ID %in% c("EM_FWOR", "EM_MAOR", "EM_PAID")) %>%
  mutate(Project_ID = str_replace(Project_ID, "_", ""))


# Select Variables of interest for model ----------------------------------

#get all variables we might want for model
locs1 <- locs %>%
  dplyr::select(Year_located, Date_located,
                Project_ID, Nest_ID, Nest_Trt,
                Transect_ID, Point_ID, Direction,
                Distance, UTM_E, UTM_N, UTM_datum_zone,
                Nest_Ht, Cavity_age, Decay_class,
                Tree_Snag_Log, Tree_sp, Tree_ht,
                DBH, Orientation, Aspect, Slope)


# Consolidate treatment IDs -----------------------------------------------

#sort facts into harvest and burn treatment categories
facts2 <- facts %>%
  mutate(category = case_when(ACTIVITY %in% c("Wildfire - Human Ignition", 
                               "Wildfire - Natural Ignition") ~ "Wildfire",
               ACTIVITY %in% c("Burning of Piled Material",
                               "Jackpot Burning - Scattered concentrations") ~ "Pile_Burn",
               ACTIVITY %in% c("Broadcast Burning - Covers a majority of the unit",
                               "Underburn - Low Intensity (Majority of Unit)",
                               "Control of Understory Vegetation- Burning") ~ "Broadcast_Burn",
               str_detect(ACTIVITY, "Thin|Cut|thinned|thinning") ~ "Harvest_Thin",
               TRUE ~ "O"))  

facts2 <- facts %>%
  mutate(category = case_when(ACTIVITY %in% c("Wildfire - Human Ignition", 
                                              "Wildfire - Natural Ignition") ~ "Wildfire",
                              ACTIVITY %in% c("Burning of Piled Material",
                                              "Jackpot Burning - Scattered concentrations") ~ "Harvest_Thin",
                              ACTIVITY %in% c("Broadcast Burning - Covers a majority of the unit",
                                              "Underburn - Low Intensity (Majority of Unit)",
                                              "Control of Understory Vegetation- Burning") ~ "Broadcast_Burn",
                              str_detect(ACTIVITY, "Thin|Cut|thinned|thinning") ~ "Harvest_Thin",
                              TRUE ~ "O"))  

#get all the unique burn treatments and last burn year
burns <- facts2 %>%
  filter(category == "Broadcast_Burn") %>%
  distinct(Nest_ID, YR_Affected, category) %>%
  group_by(Nest_ID) %>%
  summarise(n_burn = n(),
            last_burn = max(YR_Affected))

#get all the unique harvest treatments and last harvest year
harvests <- facts2 %>%
  filter(category == "Harvest_Thin") %>%
  distinct(Nest_ID, YR_Affected, category) %>%
  group_by(Nest_ID) %>%
  summarise(n_harvest = n(), #times harvested
            last_harvest = max(YR_Affected))

facts2 %>%
  distinct(Nest_ID, category) %>%
  group_by(category) %>%
  tally()
#get the categories for each nest to be "H", "B", or "HB"
categories <- facts2 %>%
  distinct(Nest_ID, category) %>%
  group_by(Nest_ID) %>%
  mutate(Trt_cat = case_when(all(category == "Broadcast_Burn") ~ "B",
                             all(category == "Harvest_Thin") ~ "H",
                             all(category == "Wildfire") ~ "W",
                             TRUE ~ "HB")) %>%
  ungroup() %>%
  distinct(Nest_ID, Trt_cat) #just get one value per nest

categories %>% 
  group_by(Trt_cat) %>%
  tally()
#combine the harvest, burn, and category DFs together so we hvae
#per treated nests 1. number of times harvest/burn/tx, 2. time since ha/bu/tx
#3. the tx category for ecah nest
all_tx <- burns %>%
  full_join(harvests, by = "Nest_ID") %>%
  #get the sum of the number of harvest and burn treatments
  mutate(n_tx = rowSums(across(c(n_harvest,n_burn)), na.rm = T)) %>%
  #get the maxiumum of the last treatment value
  mutate(last_tx =pmax(last_burn, last_harvest, na.rm = T)) %>%
  full_join(categories, by = "Nest_ID")


# Combine all tx covariates into OG DF ------------------------------------


#join that back up with the nest dataframe
locs2 <- locs1 %>%
  #join the tx variables into the nest location variable df
  full_join(all_tx, by = "Nest_ID") %>%
  ungroup() %>%
  #replace teh count data NA values with 0s
  mutate(n_burn = replace_na(n_burn, 0),
         n_harvest = replace_na(n_harvest, 0),
         n_tx = replace_na(n_tx, 0)) %>%
  #calculte a "years since treatment" by row
  rowwise() %>%
  mutate(Tm_since_b = as.character(Year_located - last_burn),
         Tm_since_h = as.character(Year_located - last_harvest),
         Tm_since_tx = as.character(Year_located - last_tx)) %>%
  mutate(Tm_since_b = replace_na(Tm_since_b, "oot"),
         Tm_since_h = replace_na(Tm_since_h, "oot"),
         Tm_since_tx = replace_na(Tm_since_tx, "oot")) %>%
  #any NA tx categories at this poitn are untreated
  mutate(Trt_cat = replace_na(Trt_cat, "U"))

locs2 %>%
  group_by(Trt_cat) %>%
  tally()

# Export ------------------------------------------------------------------

write.csv(locs2, here("data_outputs", 
                      "01_cleaning",
                      "01_Nest_trt_covs.csv"))

#END SCRIPT
