# Tree survival data prep
# November 1, 2021
# Ana Miller-ter Kuile

# prepping data for the model of tree survival

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

data <- read.csv(here('data_outputs',
                       '02_trees',
                       '02_analysis_ready',
                       'yearly_tree_data.csv'))

# Arrange by interval number for custom model -----------------------

trees2 <- data %>%
  distinct(CoreID, PlotID, Block, interval, DBH, response) %>%
  group_by(CoreID) %>%
  mutate(n.t = n()) %>%
  ungroup() %>%
  arrange(n.t, CoreID, interval)

Trees <- trees2 %>%
  distinct(CoreID, n.t)

# Count variables for indexing --------------------------------------------

# Numbers -----------------------------------------------------------------

#number with one interval only
n.trees1 <- Trees %>%
  filter(n.t == 1) %>%
  tally() %>%
  as_vector()

#total number of plants
n.trees <- length(unique(trees2$CoreID))

#How many times did each tree get measured
#(number of intervals)
n.t <- as_vector(Trees$n.t)

n.trt <- length(unique(data$Trt_cat))

n.years <- data %>%
  distinct(CoreID, interval, Year) %>%
  group_by(CoreID) %>%
  arrange(match(CoreID, Trees$CoreID)) %>%
  tally() %>%
  ungroup() %>%
  arrange(match(CoreID, Trees$CoreID)) %>%
  dplyr::select(n) %>%
  as_vector()

start.year <- data %>%
  distinct(CoreID, interval, yearID) %>%
  group_by(CoreID, interval) %>%
  filter(yearID == min(yearID)) %>%
  arrange(match(CoreID, Trees$CoreID)) %>%
  pivot_wider(names_from = 'interval',
              values_from = 'yearID') %>%
  column_to_rownames(var = "CoreID") %>%
    ungroup() %>%
  as.matrix()

end.year <- data %>%
  distinct(CoreID, interval, yearID) %>%
  group_by(CoreID, interval) %>%
  filter(yearID == max(yearID)) %>%
  arrange(match(CoreID, Trees$CoreID)) %>%
  pivot_wider(names_from = 'interval',
              values_from = 'yearID') %>%
  column_to_rownames(var = "CoreID") %>%
  ungroup() %>%
  as.matrix()

# Indexing IDs for plots, blocks, intervals -------------------------------

index_df <- data %>%
  arrange(match(CoreID, Trees$CoreID)) 

Interval.ID <- index_df %>%
  dplyr::select(CoreID, yearID, interval) %>%
  pivot_wider(names_from = 'yearID',
              values_from = 'interval') %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()


# Covariates --------------------------------------------------------------

#k's are all relative to individual tree
#Treatment ID[i, k]
TreatmentID <- index_df %>%
  mutate(Trt_cat = factor(Trt_cat, levels = c("U", "T", "TB"))) %>%
  mutate(Trt_cat = as.numeric(Trt_cat)) %>%
  dplyr::select(CoreID, yearID, Trt_cat) %>%
  pivot_wider(names_from = 'yearID',
              values_from = 'Trt_cat') %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()
  
#DBH[i, Interval.ID[i,k]]
DBH <- index_df %>%
  distinct(CoreID, DBH, interval) %>%
  mutate(DBH = scale(DBH)) %>%
  pivot_wider(names_from = 'interval',
              values_from = 'DBH') %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()

#BA[i, Interval.ID[i,k]]
BA <- index_df %>%
  distinct(CoreID, scaled_BA, interval) %>%
  pivot_wider(names_from = 'interval',
              values_from = 'scaled_BA') %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()

CanopyCover <- index_df %>%
  dplyr::select(CoreID, CC_scaled, yearID) %>%
  pivot_wider(names_from = 'yearID',
              values_from = 'CC_scaled') %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()

maxVPD <- index_df %>%
  dplyr::select(CoreID, VPD_scaled, yearID) %>%
  pivot_wider(names_from = 'yearID',
              values_from = 'VPD_scaled') %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()

minSWA <- index_df %>%
  dplyr::select(CoreID, SWA_scaled, yearID) %>%
  pivot_wider(names_from = 'yearID',
              values_from = 'SWA_scaled') %>%
  column_to_rownames(var = "CoreID") %>%
  as.matrix()

# y and y2 ----------------------------------------------------------------

y <- index_df %>%
  mutate(fate = case_when(response == "Dead" ~ 0,
                          response == "Live" ~ 1)) %>%
  group_by(CoreID) %>%
  mutate(final_fate = case_when(any(fate == 0) ~ 0,
                                TRUE ~ 1)) %>%
  distinct(CoreID, final_fate) %>%
  ungroup() %>%
  dplyr::select(final_fate) %>%
  as_vector()

y2 <- y

# Export as RDS -----------------------------------------------------------

all_data <- list(n.trees1 = n.trees1, 
                 n.trees = n.trees,
                 n.t = n.t,
                 n.trt = n.trt,
                 n.years = n.years,
                 start.year = start.year,
                 end.year = end.year,
                 Interval.ID = Interval.ID,
                 DBH = DBH,
                 BA = BA,
                 CanopyCover = CanopyCover,
                 maxVPD = maxVPD,
                 minSWA = minSWA,
                 TreatmentID = TreatmentID,
                 y = y,
                 y2 = y2)


saveRDS(all_data, here("data_outputs", 
                       "02_trees",
                       '03_JAGS_input_data',
                       "mod3_JAGS_input_data.RDS"))

saveRDS(all_data, here("monsoon", 
                       "03_trees",
                       "model3",
                       "inputs",
                       "mod3_JAGS_input_data.RDS"))

