#Compiling all covariate data for nest survival
# Ana Miller-ter Kuile
# April 4, 2023

# This script combines all the nest data with the climate data and 
#exports the final dataset taht will be used in analyses

#one will have interval-specific values, the other will just
# have final values with mean of anything that varies by interval

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "readxl", "lubridate")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

climate <- read.csv(here("data_outputs",
                         "01_whwonests",
                         "01_cleaning",
                         "05_dailyWeather_nests.csv"))

  
nests <- read.csv(here('data_outputs', 
                       "01_whwonests",
                 "01_cleaning", 
                 "04_nest_local_land_data.csv"))

# interval dataset --------------------------------------------------------

climate2 <- climate %>%
  #set climate to NA when start == end
  mutate(maxTmax_C = case_when(Julian_start != Julian_end ~ maxTmax_C,
                               Julian_start == Julian_end ~ NA_real_,
                               TRUE ~ NA_real_),
         meanTmax_C = case_when(Julian_start != Julian_end ~ meanTmax_C,
                               Julian_start == Julian_end ~ NA_real_,
                               TRUE ~ NA_real_),
         maxPpt_mm = case_when(Julian_start != Julian_end ~ maxPpt_mm,
                               Julian_start == Julian_end ~ NA_real_,
                               TRUE ~ NA_real_),
         meanPpt_mm = case_when(Julian_start != Julian_end ~ meanPpt_mm,
                               Julian_start == Julian_end ~ NA_real_,
                               TRUE ~ NA_real_)) %>%
  #get rid of starts taht match ends
  mutate(Julian_start = case_when(Julian_start == Julian_end ~ NA_integer_,
                                  TRUE ~ Julian_start)) %>%
  group_by(Nest_ID) %>%
  arrange(Julian_end) %>%
  #get intervals for each nest to combine them
  mutate(interval = 0:(n()-1)) %>%
  ungroup() %>%
  #remove rows to combine with other dataset
  dplyr::select(-Year_located, -Julian_start, -Julian_end, -UTM_datum_zone)

#remove zero-day intervals
nests_int <- nests %>%
  group_by(Nest_ID) %>%
  arrange(Julian_end) %>%
  #get same intervals for joining
  mutate(interval = 0:(n()-1)) %>%
  #join by Nest and interval
  left_join(climate2, by = c("Nest_ID", "interval")) %>%
  #select important covariates only
  dplyr::select(Nest_ID, Fate_cat, Visit_date,
                Year_located, Project_ID, Transect_ID2, 
                Julian_start, Julian_end, interval,
                Stage, 
                Trt_cat,  Nest_Ht, cosOrientation, Init_day, Tree_sp,
                Trees_2550, Trees_50, pPIPO,
                a1000_areacv2, a1000_contag, a1000_np1,
                a1000_Ha, a1000_RxBu,
                maxTmax_C, meanTmax_C, maxPpt_mm, meanPpt_mm
                )

#some nests have zero-day intervals, so we would like to get rid of those

first_obs <- nests_int %>%
  group_by(Nest_ID) %>%
  #find the minimum value per nest
  mutate(min_end = min(Julian_end, na.rm = T)) %>%
  #set any matching start-end that is the minimum (first) for that nest to NA
  mutate(Julian_start = case_when(Julian_start == Julian_end & Julian_start == min_end ~ NA_integer_,
                                  TRUE ~ Julian_start)) %>%
  ungroup() %>%
  dplyr::select(-min_end) 
  
#find those nests that now still have zero-day intervals
zero_ints <- first_obs %>%
    filter(Julian_start == Julian_end)

#remove zero-day intervals from dataset
nests_int2 <- first_obs %>%
  anti_join(zero_ints) %>%
  group_by(Nest_ID) %>%
  mutate(prevStage = lag(Stage)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(Age = Julian_end - Init_day) %>%
  filter(!prevStage %in% c("E", "F")) %>%
  #remove any that didn't have a previous stage (time zero obs)
  filter(!is.na(prevStage)) %>%
  #make NA julian_start == initiation day
  mutate(Julian_start = case_when(is.na(Julian_start) ~ Init_day,
                                  TRUE ~ Julian_start))

write.csv(nests_int2, here("data_outputs",
                           "01_whwonests",
               "02_analysis_ready",
               "01_whwonests",
               "interval_models_nest_data.csv"))

# Full-survey model data --------------------------------------------------

#Get total exposure time by getting the initition day (when known) and the final
# survey day. If inititaion date is unknown, set start to be the start of the first
#survey interval

t <- nests_int2 %>%
  group_by(Nest_ID) %>%
  mutate(min = Init_day,
         max = max(Julian_end, na.rm = T)) %>%
  mutate(min = case_when(is.na(min) ~ min(Julian_start),
                         TRUE ~ min)) %>%  
  ungroup() %>%
  rowwise() %>%
  #get total days to be max - min
  mutate(exposure = max - min) %>%
  #select these to join below
  distinct(Nest_ID, exposure)

stage <- nests_int2 %>%
  group_by(Nest_ID) %>%
  filter(interval == max(interval, na.rm = T)) %>%
  dplyr::select(Nest_ID, prevStage)

#get by nest values
nests_full <- nests_int2 %>%
  group_by(Nest_ID, Fate_cat, Year_located, Project_ID, Transect_ID2,
           Trt_cat, Nest_Ht, cosOrientation, Init_day, Tree_sp,
           Trees_2550, Trees_50, pPIPO, 
           a1000_areacv2, a1000_contag, a1000_np1, a1000_Ha, a1000_RxBu) %>%
  summarise(maxTmax_C = max(maxTmax_C, na.rm = T),
         meanTmax_C = mean(meanTmax_C, na.rm = T),
         maxPpt_mm = max(maxPpt_mm, na.rm = T),
         meanPpt_mm = mean(meanPpt_mm, na.rm = T)) %>%
  ungroup() %>%
  left_join(t, by = "Nest_ID") %>%
  left_join(stage, by = "Nest_ID")

write.csv(nests_full, here("data_outputs",
                           "01_whwonests",
                           "02_analysis_ready",
                           "01_whwonests",
                           "full_survey_model_nest_data.csv"))



