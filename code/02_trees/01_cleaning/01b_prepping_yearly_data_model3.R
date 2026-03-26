
package.list <- c("here", "tidyverse" , "lubridate")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

trees <- read.csv(here("data_raw",
                       "02_trees",
                       "treeMortality.csv"))

climate <- read.csv(here('data_raw',
                         '02_trees',
                         'climateSummaries_seasonal.csv'))

canopy <- read.csv(here('data_raw',
                        '02_trees',
                        'annualCanopyCover.csv'))



# Subset site FV from data ------------------------------------------------

#just chooosing fort valley trees
trees2 <- trees %>%
  filter(Site %in% c("GV", "AS"))

# Get yearly climate per plot ---------------------------------------------

#weight things by season length (how many months)
#dry summer - 2
#fall -2
#monsoon summer -3
#spring - 2
#winter - 3

#take min or max of the annual *means*

climate2 <- climate %>%
  rowwise() %>%
  mutate(meanVPD = weighted.mean(x = c(meanVPD_hPa_Dry.Summer, meanVPD_hPa_Fall,
                                       meanVPD_hPa_Monsoon.Summer, meanVPD_hPa_Spring,
                                       meanVPD_hPa_Winter),
                                 w = c(2, 2, 3, 2, 3), na.rm = T),
         meanSWA = weighted.mean(x = c(SWA_mm_000to150_cm_Dry.Summer, SWA_mm_000to150_cm_Fall,
                                       SWA_mm_000to150_cm_Monsoon.Summer, SWA_mm_000to150_cm_Spring,
                                       SWA_mm_000to150_cm_Winter), 
                                 w =c(2, 2, 3, 2, 3), na.rm = T)) %>%
  dplyr::select(PlotID, Thinned, Burned, Year, meanVPD,
                meanSWA) %>%
  ungroup() 

climate3 <- climate2 %>%
  filter(str_detect(PlotID, "GV|AS"))

canopy2 <- canopy %>%
  filter(str_detect(PlotID, "GV|AS"))


# "Stretch" trees to all yars ---------------------------------------------

years_levels <- as.factor(1997:2019)

trees3 <- trees2 %>%
  group_by(CoreID) %>%
  arrange(Year) %>%
  mutate(interval = 1:n()) %>%
  ungroup()

trees4 <- trees3 %>%
  distinct(CoreID, PlotID, Block, priorDBH, priorBA, 
           Year, priorYear, interval, response) %>%
  #subtracting a year from previous survey since it 
  #survived through that year and the survey year is 
  #still ongoing during survey, so include in next interval
  mutate(Year = case_when(Year < 2019 ~ Year - 1,
                               TRUE ~ Year)) %>%
  group_by(CoreID) %>%
  # mutate(minYear = min(priorYear),
  #        maxYear = max(Year)) %>%
  mutate(year = map2(priorYear, Year,
                      ~seq(.x, .y, by = 1))) %>%
  unnest(year) %>%
  ungroup() %>%
  distinct(CoreID, PlotID, Block, priorDBH, priorBA,
           year, interval, response) %>%
  group_by(CoreID) %>%
  mutate(yearID = 1:n()) %>%
  ungroup() %>%
  #rename priorDBH 
  rename("DBH" = 'priorDBH',
         "Year" = 'year')

# All yearly covariates ---------------------------------------------------

#need: treatment, DBH, basal area,
#canopy cover, vpd and swa for each year

#treatment is at the plot level and 
# at the interval level
# b1[TreatmentID[Plot.ID[i],Interval.ID[k]]] +
# #DBH is basd on the tree and that interval
# b[2]*DBH[i,Interval.ID[k]] +
# b[3]*DBH[i,Interval.ID[k]]^2 +
# #basal area is dependent on the plot
# #and interval
# b[4]*BA[Plot.ID[i],Interval.ID[k]] +
# #canopy cover is based on the plot
# #but we have yearly info
# b[5]*CanopyCover[Plot.ID[i],k] +
# #vpd is at the block level and yearly
# b[6]*maxVPD[Block.ID[i],k] +
# #soil water is at the plot level and yearly
# b[7]*minSWA[Plot.ID[i],k] 

#things at the block level:
#vpd

#things at the plot level:
#treatment
#basal area
#canopy cover
#swa?

#things at the tree level:
#DBH

#in the model, can have those
#that are at block and plot level vary at that 
#level instead of at the tree

#Tree ID
#Plot ID
#Block ID
#Interval ID
#Interval number for each tree
#years corresponding to their interval ID for each tree/plot/block

# Generate main DF --------------------------------------------------------

#want a yearly dataframe with tree ID,
#covariates
#plot, block, intervals,
#year, 

#survey years are 1997-2019
#intervals are 1997-2004, 
#2004-2008
#2008-2019

climate4 <- climate3 %>%
  filter(Year %in% c(1997:2019)) %>%
  mutate(Trt_cat = case_when((Thinned == 0 & Burned == 0) ~ "U",
                             (Thinned == 0 & Burned == 1) ~ "B",
                             (Thinned == 1 & Burned == 0) ~ "T",
                             (Thinned == 1 & Burned == 1) ~ "TB"))

cc <- canopy2 %>%
   filter(Year %in% 1997:2019) %>%
  dplyr::select(PlotID, Year, CanopyCover)

env_covs <- cc %>%
  left_join(climate4, by = c("PlotID", "Year")) %>%
  dplyr::select(PlotID, Year,  
                CanopyCover, meanVPD, meanSWA,
                Trt_cat)

tree_plot_ids <- trees3 %>%
  distinct(Year, PlotID) 

#seems like there are some plots that just don't have
#tracked trees??
env_covs2 <- env_covs %>%
  filter(PlotID %in% tree_plot_ids$PlotID) %>%
  mutate(CC_scaled = scale(CanopyCover),
         VPD_scaled = scale(meanVPD),
         SWA_scaled = scale(meanSWA))

ba <- trees4 %>%
  distinct(PlotID, interval, priorBA) %>%
  mutate(scaled_BA = scale(priorBA))

all_dat <- trees4 %>%
  left_join(env_covs2, by = c("PlotID", "Year")) %>%
  left_join(ba, by = c("PlotID", "interval", "priorBA"))

write.csv(all_dat, 
          here("data_outputs",
               '02_trees',
               "02_analysis_ready",
               "yearly_tree_data.csv"))

