#need:
#summarised TOTAL average values for all covariates for model 1
#time intervals, t, in each dataset



# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse" , "lubridate")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

trees <- read.csv(here("data_raw",
                       "03_trees",
                       "treeMortality.csv"))

climate <- read.csv(here('data_raw',
                         '03_trees',
                         'climateSummaries_seasonal.csv'))

canopy <- read.csv(here('data_raw',
                        '03_trees',
                        'annualCanopyCover.csv'))



# Subset site FV from data ------------------------------------------------

#just chooosing fort valley trees
trees2 <- trees %>%
  filter(Site == "FV")


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

# Get interval data -------------------------------------------------------

#to get the moving window of means per interval (will also repeat
# for total non-interval data)
subset_fun <- function(df, plot, start_date, end_date, fun, variable){
  
  dat <- df %>%
    filter(PlotID == plot) %>%
    filter(Year >= start_date & Year < end_date) %>%
    summarise(var = fun({.data[[variable]]})) %>%
    as_vector()
  
  return(dat)
  
  
}

#testing
subset_fun(df = climate2,
           plot = "FV_1_1_1",
           start_date = trees2$priorYear[1],
           end_date = trees2$Year[1],
           fun = mean,
           variable = "meanVPD")
  

# Interval data for canopy and climate ------------------------------------


trees2$CanopyCover <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$CanopyCover[i] <- subset_fun(df = canopy,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      fun = mean,
                                      variable = "CanopyCover")
}


trees2$maxVPD <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$maxVPD[i] <- subset_fun(df = climate2,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      fun = max,
                                      variable = "meanVPD")
}

trees2$minSWA <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$minSWA[i] <- subset_fun(df = climate2,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      fun = min,
                                      variable = "meanSWA")
}

write.csv(trees2, here("data_outputs",
                       "03_trees",
                       "02_analysis_ready",
                       "interval_tree_data.csv"))

# Total interval data -----------------------------------------------------

trees3 <- trees2 %>%
  group_by(Site, Block, Trt, Plot, Tree_Number, CoreID, PlotID) %>%
  mutate(response = as.factor(response)) %>%
  mutate(response = case_when(any(response == "Dead") ~ "Dead",
                                  TRUE ~ "Live")) %>%
  ungroup() %>%
  group_by(Site, Block, Trt, Plot, Tree_Number, CoreID, PlotID, response) %>%
  summarise(start_year = min(priorYear),
            end_year = max(Year),
            BA = mean(priorBA, na.rm = T),
            DBH = mean(priorDBH, na.rm = T)) %>%
  ungroup()



trees3$CanopyCover <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$CanopyCover[i] <- subset_fun(df = canopy,
                                      plot = trees3$PlotID[i],
                                      start_date = trees3$start_year[i],
                                      end_date = trees3$end_year[i],
                                      fun = mean,
                                      variable = "CanopyCover")
}

trees3$maxVPD <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$maxVPD[i] <- subset_fun(df = climate2,
                                             plot = trees3$PlotID[i],
                                             start_date = trees3$start_year[i],
                                             end_date = trees3$end_year[i],
                                             fun = max,
                                             variable = "meanVPD")
}

trees3$minSWA <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$minSWA[i] <- subset_fun(df = climate2,
                                         plot = trees3$PlotID[i],
                                         start_date = trees3$start_year[i],
                                         end_date = trees3$end_year[i],
                                         fun = min,
                                         variable = "meanSWA")
}


write.csv(trees3, here("data_outputs",
                       "03_trees",
                       "02_analysis_ready",
                       "total_tree_data.csv"))

