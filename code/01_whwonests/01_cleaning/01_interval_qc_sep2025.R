
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

visits <- read_xlsx(here('data_raw',
                        '01_whwonests',
                        'bird_data',
                        'Birds02_nest_visits.xlsx')) %>%
  filter(str_detect(Nest_ID, "EMFWOR|EMPAID|EMMAOR")) %>%
  dplyr::select(Nest_ID, Visit_date,
                No_eggs, No_young,
                Stage)

fates <- read_xlsx(here('data_raw',
                        '01_whwonests',
                        'bird_data',
                        'Birds03_nest_fates.xlsx')) %>%
  #cateogrize into success/failure, and unknown
  mutate(Fate_cat = case_when(Fate == 1 ~ "success",
                              Fate == 11 ~ "unknown",
                              TRUE ~ "failure")) %>%
  #remove unknown for future analyses
  filter(Fate_cat != "unknown") %>%
  filter(str_detect(Nest_ID, "EMFWOR|EMPAID|EMMAOR")) %>%
  dplyr::select(Nest_ID, Initiation_date,
                NoFL_uncert, NoFL_cert,
                Fate_cat)


# Join datasets to keep known fate ----------------------------------------

visits2 <- visits %>%
  filter(Nest_ID %in% fates$Nest_ID) %>%
  left_join(fates, by = "Nest_ID") %>%
  #remove excavating
  filter(Stage != "E")


# Get interval start-stops ------------------------------------------------

intervals <- visits2 %>%
  group_by(Nest_ID) %>%
  arrange(Nest_ID, Visit_date, No_eggs) %>%
  distinct(Nest_ID, Visit_date, No_eggs,        
           No_young, Stage, Initiation_date,
           NoFL_uncert, NoFL_cert, Fate_cat) %>%
  mutate(prev_survey = lag(Visit_date)) %>%
  ungroup()

# QC repeat visits --------------------------------------------------------

#some nests were visited twice in one day

intervalqc <- intervals %>%
  group_by(Nest_ID) %>%
  filter(as_date(Visit_date) == as_date(prev_survey)) %>%
  ungroup()

intervalsdup <- intervalqc %>%
  distinct(Nest_ID)

intervalsdup_df <- intervals %>%
  filter(Nest_ID %in% intervalsdup$Nest_ID)

intervalsdup <- intervalsdup_df %>%
  group_by(Nest_ID) %>%
  filter(!(as_date(Visit_date) == as_date(prev_survey) & No_eggs == 999)) %>%
  #still funky
  #EMMAOR_NPAWT1-03A_2021
  filter(!((Nest_ID == "EMMAOR_NPAWT1-03A_2021") & 
             (Visit_date == "2021-07-09") &
             (Stage == "N")))
             
# Update full dataset -----------------------------------------------------

intervals2 <- intervals %>%
  group_by(Nest_ID) %>%
  filter(!(as_date(Visit_date) == as_date(prev_survey) & No_eggs == 999)) %>%
  #still funky
  #EMMAOR_NPAWT1-03A_2021
  filter(!((Nest_ID == "EMMAOR_NPAWT1-03A_2021") & 
             (Visit_date == "2021-07-09") &
             (Stage == "N"))) %>%
  ungroup() %>%
  dplyr::select(-prev_survey) %>%
  group_by(Nest_ID) %>%
  arrange(Nest_ID, Visit_date) %>%
  mutate(interval = 0:(n()-1)) %>%
  mutate(Prev_date = lag(Visit_date)) %>%
  ungroup() %>%
  dplyr::select(Nest_ID, Visit_date,
                Prev_date, interval,
                Fate_cat,
                Initiation_date,
                Stage) %>%
  filter(interval > 0)


# Export ------------------------------------------------------------------

write.csv(intervals2,
        here("data_outputs",
             "01_whwonests",
             "01_cleaning",
             "01_surveys_fates_qc.csv"),
        row.names = F)  

