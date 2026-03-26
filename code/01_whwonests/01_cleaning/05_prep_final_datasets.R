#this summarizes the total and the interval models
#the daily custom model can be prepped from these 
#and other existing datasets


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "readxl")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#custom functions 
source(here("code", 
            "00_functions", 
            "tidy_functions.R"))

# Load datasets -----------------------------------------------------------

visits <- read.csv(here("data_outputs",
                        "01_whwonests",
                        "01_cleaning",
                        "02_surveys_fates_locations.csv"))

climate <- read.csv(here("data_outputs",
                         "01_whwonests",
                         "01_cleaning",
                         "03_dailyWeather_nests.csv"))

trees <- read.csv(here("data_outputs",
                       "01_whwonests",
                       "01_cleaning",
                       "04_nest_tree_densities.csv"))

# Total exposure model ----------------------------------------------------

total_nests <- visits %>%
  group_by(Nest_ID) %>%
  filter(interval == max(interval)) %>%
  dplyr::select(Nest_ID, Fate_cat,
                Initiation_date, Year_located,
                Project_ID, Nest_Ht, Orientation) %>%
  left_join(trees, by = c("Nest_ID" = "Measurement_ID")) %>%
  dplyr::select(Nest_ID, Fate_cat, 
                Initiation_date, Year_located,
                Project_ID, Nest_Ht,
                Orientation, Trees_2550,
                Trees_50)

total_climate <- climate %>%
  group_by(Nest_ID) %>%
  reframe(tmax = mean(tmax, na.rm = T),
          ppt = mean(ppt, na.rm = T)) %>%
  ungroup()

exposure <- visits %>%
  group_by(Nest_ID) %>%
  reframe(first = min(Prev_date),
          last = max(Visit_date)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(t = as_date(last) - as_date(first)) %>%
  dplyr::select(-last, -first) %>%
  ungroup()
  

total_nests2 <- total_nests %>%
  left_join(total_climate, by = "Nest_ID") %>%
  left_join(exposure, by = "Nest_ID")

# Interval nests ----------------------------------------------------------

nest_intervals <- visits %>%
  dplyr::select(Nest_ID, Visit_date, Prev_date,
                interval) %>%
  rowwise() %>%
  mutate(Prev_date = case_when(interval > 1 ~ as_date(Prev_date)+1,
                               TRUE ~ as_date(Prev_date))) %>%
  mutate(Visit_date = as_date(Visit_date)) %>%
  pivot_longer(Visit_date:Prev_date,
               names_to = "period",
               values_to = "Date") %>%
  dplyr::select(-period) %>%
  distinct(Nest_ID, interval, Date)

climate_intervals <- climate %>%
  mutate(Date = as_date(Date)) %>%
  left_join(nest_intervals, by = c("Nest_ID", "Date")) %>%
  fill(interval, .direction = 'down')

climate_sum_intervals <- climate_intervals %>%
  group_by(Nest_ID, interval) %>%
  reframe(tmax = mean(tmax, na.rm = T),
          ppt = mean(ppt, na.rm = T))

interval_data <- visits %>%
  left_join(climate_sum_intervals, by = c("Nest_ID", "interval")) %>%
  left_join(trees, by = c("Nest_ID" = "Measurement_ID")) %>%
  dplyr::select(Nest_ID, Fate_cat, Visit_date, Prev_date, interval,
                Initiation_date, Year_located,
                Project_ID, Nest_Ht,
                Orientation, Trees_2550,
                Trees_50, tmax, ppt)


# Export datasets ---------------------------------------------------------

write.csv(total_nests2, 
          here("data_outputs",
               "01_whwonests",
               "02_analysis_ready",
               'full_survey_model_nest_data.csv'),
          row.names = F)

write.csv(climate_intervals, 
          here("data_outputs",
               "01_whwonests",
               "02_analysis_ready",
               'daily_climate_data_withintervalID.csv'),
          row.names = F)

write.csv(interval_data, 
          here("data_outputs",
               "01_whwonests",
               "02_analysis_ready",
               'interval_models_nest_data.csv'),
          row.names = F)
