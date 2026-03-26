
# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse", 
                  'lubridate', "readxl")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

#what we need:
#intervals and fates data
#location data that has info on nests

visits <- read.csv(here("data_outputs",
                        "01_whwonests",
                        "01_cleaning",
                        "01_surveys_fates_qc.csv"))

locations <- read_xlsx(here("data_raw",
                           '01_whwonests',
                           'bird_data',
                           'Birds01_nest_locations.xlsx')) %>%
  dplyr::select(Year_located, Project_ID,
                Nest_ID, UTM_E, UTM_N,
                UTM_datum_zone,
                Nest_Ht, Orientation) %>%
  filter(Project_ID %in% c("EM_FWOR",
                           "EM_MAOR", 
                           "EM_PAID")) %>%
  filter(Nest_ID %in% visits$Nest_ID)


# Join datasets -----------------------------------------------------------

visits2 <- visits %>%
  left_join(locations, by = "Nest_ID") %>%
  dplyr::select(-X)


# Export ------------------------------------------------------------------

write.csv(visits2,
          here("data_outputs",
               "01_whwonests",
               "01_cleaning",
               "02_surveys_fates_locations.csv"),
          row.names = F)
