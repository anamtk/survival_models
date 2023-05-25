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

# Get interval data -------------------------------------------------------

#to get the moving window of means per interval (will also repeat
# for total non-interval data)
subset_fun <- function(df, plot, start_date, end_date, variable){
  
  dat <- df %>%
    filter(PlotID == plot) %>%
    filter(Year >= start_date & Year < end_date)
  
  mean(dat[,variable], na.rm = TRUE)
  
}

#testing
subset_fun(df = canopy,
           plot = "FV_1_1_1",
           start_date = trees2$priorYear[1],
           end_date = trees2$Year[1],
           variable = "CanopyCover")
  

# Interval data for canopy and climate ------------------------------------


trees2$CanopyCover <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$CanopyCover[i] <- subset_fun(df = canopy,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "CanopyCover")
}

trees2$meanVPD_Dry.Summer <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$meanVPD_Dry.Summer[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "meanVPD_hPa_Dry.Summer")
}

trees2$meanVPD_Fall <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$meanVPD_Fall[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "meanVPD_hPa_Fall")
}

trees2$meanVPD_Monsoon.Summer <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$meanVPD_Monsoon.Summer[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "meanVPD_hPa_Monsoon.Summer")
}

trees2$meanVPD_Spring <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$meanVPD_Spring[i] <- subset_fun(df = climate,
                                                 plot = trees2$PlotID[i],
                                                 start_date = trees2$priorYear[i],
                                                 end_date = trees2$Year[i],
                                                 variable = "meanVPD_hPa_Spring")
}

trees2$meanVPD_Winter <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$meanVPD_Winter[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "meanVPD_hPa_Winter")
}

trees2$SWA_Dry.Summer <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$SWA_Dry.Summer[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "SWA_mm_000to150_cm_Dry.Summer")
}

trees2$SWA_Fall <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$SWA_Fall[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "SWA_mm_000to150_cm_Fall")
}

trees2$SWA_Monsoon.Summer <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$SWA_Monsoon.Summer[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "SWA_mm_000to150_cm_Monsoon.Summer")
}

trees2$SWA_Spring <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$SWA_Spring[i] <- subset_fun(df = climate,
                                         plot = trees2$PlotID[i],
                                         start_date = trees2$priorYear[i],
                                         end_date = trees2$Year[i],
                                         variable = "SWA_mm_000to150_cm_Spring")
}

trees2$SWA_Winter <- rep(NA, nrow(trees2))

for(i in 1:nrow(trees2)){
  trees2$SWA_Winter[i] <- subset_fun(df = climate,
                                      plot = trees2$PlotID[i],
                                      start_date = trees2$priorYear[i],
                                      end_date = trees2$Year[i],
                                      variable = "SWA_mm_000to150_cm_Winter")
}

write.csv(trees2, here("data_outputs",
                       "03_trees",
                       "interval_tree_data.csv"))

# Total interval data -----------------------------------------------------

trees3 <- trees2 %>%
  group_by(Site, Block, Trt, Plot, Tree_Number, CoreID, PlotID) %>%
  summarise(start_year = min(Year),
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
                                      variable = "CanopyCover")
}

trees3$meanVPD_Dry.Summer <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$meanVPD_Dry.Summer[i] <- subset_fun(df = climate,
                                             plot = trees3$PlotID[i],
                                             start_date = trees3$start_year[i],
                                             end_date = trees3$end_year[i],
                                             variable = "meanVPD_hPa_Dry.Summer")
}

trees3$meanVPD_Fall <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$meanVPD_Fall[i] <- subset_fun(df = climate,
                                       plot = trees3$PlotID[i],
                                       start_date = trees3$start_year[i],
                                       end_date = trees3$end_year[i],
                                       variable = "meanVPD_hPa_Fall")
}

trees3$meanVPD_Monsoon.Summer <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$meanVPD_Monsoon.Summer[i] <- subset_fun(df = climate,
                                                 plot = trees3$PlotID[i],
                                                 start_date = trees3$start_year[i],
                                                 end_date = trees3$end_year[i],
                                                 variable = "meanVPD_hPa_Monsoon.Summer")
}

trees3$meanVPD_Spring <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$meanVPD_Spring[i] <- subset_fun(df = climate,
                                         plot = trees3$PlotID[i],
                                         start_date = trees3$start_year[i],
                                         end_date = trees3$end_year[i],
                                         variable = "meanVPD_hPa_Spring")
}

trees3$meanVPD_Winter <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$meanVPD_Winter[i] <- subset_fun(df = climate,
                                         plot = trees3$PlotID[i],
                                         start_date = trees3$start_year[i],
                                         end_date = trees3$end_year[i],
                                         variable = "meanVPD_hPa_Winter")
}

trees3$SWA_Dry.Summer <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$SWA_Dry.Summer[i] <- subset_fun(df = climate,
                                         plot = trees3$PlotID[i],
                                         start_date = trees3$start_year[i],
                                         end_date = trees3$end_year[i],
                                         variable = "SWA_mm_000to150_cm_Dry.Summer")
}

trees3$SWA_Fall <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$SWA_Fall[i] <- subset_fun(df = climate,
                                   plot = trees3$PlotID[i],
                                   start_date = trees3$start_year[i],
                                   end_date = trees3$end_year[i],
                                   variable = "SWA_mm_000to150_cm_Fall")
}

trees3$SWA_Monsoon.Summer <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$SWA_Monsoon.Summer[i] <- subset_fun(df = climate,
                                             plot = trees3$PlotID[i],
                                             start_date = trees3$start_year[i],
                                             end_date = trees3$end_year[i],
                                             variable = "SWA_mm_000to150_cm_Monsoon.Summer")
}

trees3$SWA_Spring <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$SWA_Spring[i] <- subset_fun(df = climate,
                                     plot = trees3$PlotID[i],
                                     start_date = trees3$start_year[i],
                                     end_date = trees3$end_year[i],
                                     variable = "SWA_mm_000to150_cm_Spring")
}

trees3$SWA_Winter <- rep(NA, nrow(trees3))

for(i in 1:nrow(trees3)){
  trees3$SWA_Winter[i] <- subset_fun(df = climate,
                                     plot = trees3$PlotID[i],
                                     start_date = trees3$start_year[i],
                                     end_date = trees3$end_year[i],
                                     variable = "SWA_mm_000to150_cm_Winter")
}


write.csv(trees3, here("data_outputs",
                       "03_trees",
                       "total_tree_data.csv"))

