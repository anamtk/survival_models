#Compiling all covariate data for kelp survival
# Ana Miller-ter Kuile
# April 27, 2023

# This script takes interval kelp data and summarises to 
# the total survey period

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse" , "lubridate")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

kelp <- read.csv(here("data_raw",
                      "02_kelp",
                      "Kelp_static_visits_sst_waves_substr_slope.csv"))



# Remove two plants that were only surveyed once --------------------------

kelp_int <- kelp %>%
  filter(!Card_number %in% c(2407, 2505))

write.csv(kelp_int, here("data_outputs",
                           "02_kelp",
                         "02_analysis_ready",
                           "interval_kelp_data.csv"))

# Create total survey dataset ---------------------------------------------

kelp_total <- kelp_int %>%
  group_by(Card_number, 
           Site_transect, Site, Transect, 
           Depth_ft, Max_Diam_1_cm,
           Substrate_Code) %>%
  summarise(Stipe_num = mean(Stipe_num, na.rm = T),
            sst_mean = mean(sst_mean, na.rm = T),
            wave_p_mean = mean(wave_p_mean, na.rm = T)) %>%
  ungroup()

t <- kelp_int %>%
  group_by(Card_number) %>%
  summarise(max = max(as.Date(Visit_Date), na.rm = T),
         min = min(as.Date(Visit_Date), na.rm = T)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(t = max - min) %>%
  dplyr::select(-max, -min)

presence <- kelp_int %>%
  dplyr::select(Card_number, Visit_interval, Presence) %>%
  group_by(Card_number) %>%
  filter(Visit_interval == max(Visit_interval, na.rm = T)) %>%
  ungroup() %>%
  dplyr::select(-Visit_interval)

kelp_total <- kelp_total %>%
  left_join(t, by = "Card_number") %>%
  left_join(presence, by = "Card_number")
  
write.csv(kelp_total, here("data_outputs",
                           "02_kelp",
                           "02_analysis_ready",
                           "total_intervals_kelp_data.csv"))
