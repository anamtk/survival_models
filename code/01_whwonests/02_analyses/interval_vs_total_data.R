#Explornig variation in nest climate data
# Ana Miller-ter Kuile
#April 27, 2023

#this script explores variation in the climate data for the 
#survival mdoel applied to the nest data - will repeat this exploration
# for other datasets as well

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "patchwork", "glmmTMB")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

# Load data ---------------------------------------------------------------

nests_int2 <- read.csv(here("data_outputs",
                            "02_analysis_ready",
                            "01_whwonests",
                            "interval_models_nest_data.csv"))


nests_full <- read.csv(here("data_outputs",
     "02_analysis_ready",
     "01_whwonests",
     "full_survey_model_nest_data.csv"))
  

# Histogram distributions  ------------------------------------------------

#full dataset
(a <- ggplot(nests_int2, aes(x = meanTmax_C)) +
  geom_histogram() +
  xlim(10, 39))

(b <- ggplot(nests_int2, aes(x = meanPpt_mm)) +
    geom_histogram() +
    xlim(-0.5, 9.6))

#summarised
c <- ggplot(nests_full, aes(x = meanTmax_C)) +
  geom_histogram() +
  xlim(10, 39)

d <- ggplot(nests_full, aes(x = meanPpt_mm)) +
  geom_histogram() +
  xlim(-0.5, 9.6)

#plot
(a +b)/(c+d)


# Correlation b/w interval and total data ---------------------------------

tnest <- nests_int2 %>%
  dplyr::select(Nest_ID, meanTmax_C, meanPpt_mm)

t2 <- nests_full %>%
  dplyr::select(Nest_ID, meanTmax_C, meanPpt_mm) %>%
  rename('totTmax' = 'meanTmax_C',
         'totPPT' = 'meanPpt_mm')

tnest <- tnest %>%
  left_join(t2, by = "Nest_ID")

#temp correlation
tmean <- ggplot(tnest, aes(x = totTmax, y = meanTmax_C)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_smooth(method = "lm", se =F) +
  labs(x = "Total interval mean temperature",
       y = "Interval mean temperature")

#ppt correlation
ppt <- ggplot(tnest, aes(x = totPPT, y = meanPpt_mm)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Total interval mean precipitation",
       y = "Interval mean precipitation")

#plot it
tmean + ppt

# Temperature and PPT through the season ----------------------------------

m1 <- glmmTMB(meanTmax_C ~ Julian_end + (1|Nest_ID),
              data = nests_int2)

summary(m1)
ggplot(nests_int2, aes(x = Julian_end, y = meanTmax_C, group = Nest_ID)) +
  geom_point() +
  geom_abline(slope = 0.21, intercept = -13.48)

m2 <- glmmTMB(meanPpt_mm ~ Julian_end + (1|Nest_ID),
              data = nests_int2)

summary(m2)

ggplot(nests_int2, aes(x = Julian_end, y = meanPpt_mm)) +
  geom_point() +
  geom_abline(slope = -0.02, intercept = 3.99)
