# Tree survival covariate correlation/selection
# November 1, 2021
# Ana Miller-ter Kuile

# prepping data for the model of nest survival

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "MuMIn",
                  "GGally")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

trees <- read.csv(here("data_outputs",
                       "03_trees",
                       "02_analysis_ready",
                       "interval_tree_data.csv"))


# Correlation check -------------------------------------------------------

vpd_check <- trees %>%
  dplyr::select("meanVPD_Dry.Summer", "meanVPD_Fall",          
                "meanVPD_Monsoon.Summer", "meanVPD_Spring",        
                "meanVPD_Winter")

ggpairs(vpd_check)

#fall - dry summer
#spring - dry summer
#winter - spring
swa_check <- trees %>%
  dplyr::select("SWA_Dry.Summer",        
                "SWA_Fall", "SWA_Monsoon.Summer",    
                "SWA_Spring", "SWA_Winter")

ggpairs(swa_check)

#monsoon summer - fall
#winter - spring


# VPD ---------------------------------------------------------------------
trees2 <- trees %>%
  mutate(y = case_when(response == "Live" ~ 1,
                       response == "Dead" ~ 0))

#fall - dry summer
#spring - dry summer
#winter - spring

m1 <- glm(y ~ meanVPD_Dry.Summer,
         data = trees2,
         family = "binomial")

m2 <- update(m1,  ~ -meanVPD_Dry.Summer + meanVPD_Fall)

m3 <- update(m1,  ~ -meanVPD_Dry.Summer + meanVPD_Spring)

m4 <- update(m1,  ~ -meanVPD_Dry.Summer + meanVPD_Winter)

AICc(m1, m2, m3, m4)
#fall
#dry.summer
#spring
#winter

#since fall corr with dry.summer, choose
#FALL
#since spring is better than winter
#SPRING

# SWA ---------------------------------------------------------------------

#monsoon summer - fall
#winter - spring

m5 <- glm(y ~ SWA_Monsoon.Summer,
          data = trees2,
          family = "binomial")

m6 <- update(m5,  ~ -SWA_Monsoon.Summer + SWA_Fall)

m7 <- update(m5,  ~ -SWA_Monsoon.Summer + SWA_Winter)

m8 <- update(m5,  ~ -SWA_Monsoon.Summer + SWA_Spring)

AICc(m5, m6, m7, m8)

#fall
#monsoon
#winter
#spring

#Choose:
#FALL & WINTER



