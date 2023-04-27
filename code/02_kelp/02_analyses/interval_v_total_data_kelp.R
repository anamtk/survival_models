#Explornig variation in kelp interval data
# Ana Miller-ter Kuile
#April 27, 2023

#this script explores variation in the interval covariates for the 
#survival mdoel applied to the kelp data - will repeat this exploration
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

kelp <- read.csv(here("data_outputs",
                      "02_kelp",
                      "Kelp_static_visits_sst_waves_substr_slope.csv"))


# Create total survey dataset ---------------------------------------------

kelp_total <- kelp %>%
  group_by(Card_number) %>%
  summarise(meanstipe = mean(prev_Stipes, na.rm = T),
            meansst = mean(sst_mean, na.rm = T),
            meanwave = mean(wave_p_mean, na.rm = T)) %>%
  ungroup()


# Histograms of distributions ---------------------------------------------

a <- ggplot(kelp, aes(x = prev_Stipes)) +
  geom_histogram() +
  xlim(-5, 425)

b <- ggplot(kelp, aes(x = sst_mean)) +
  geom_histogram() +
  xlim(12, 22)

c <- ggplot(kelp, aes(x = wave_p_mean)) +
  geom_histogram() +
  xlim(0, 6)

d <- ggplot(kelp_total, aes(x = meanstipe)) +
  geom_histogram() +
  xlim(-5, 425)

e <- ggplot(kelp_total, aes(x = meansst)) +
  geom_histogram() +
  xlim(12, 22)

f <- ggplot(kelp_total, aes(x = meanwave)) +
  geom_histogram() +
  xlim(0, 6)

(a + b+ c)/(d + e+ f)


# Correlation between interval and total ----------------------------------

kelp2 <- kelp %>%
  dplyr::select(Card_number, prev_Stipes, 
                sst_mean, wave_p_mean) %>%
  left_join(kelp_total, by = "Card_number")

ggplot(kelp2, aes(x = meanstipe, y = prev_Stipes)) +
  geom_abline(slope = 1) +
  geom_point() +
  geom_smooth(method = "lm", se =F) 

ggplot(kelp2, aes(x = meansst, y = sst_mean)) +
  geom_abline(slope = 1) +
  geom_point() +
  geom_smooth(method = "lm", se =F) 

ggplot(kelp2, aes(x = meanwave, y = wave_p_mean)) +
  geom_abline(slope = 1) +
  geom_point() +
  geom_smooth(method = "lm", se =F) 


# Seasonal effects --------------------------------------------------------

#do kelp get bigger as they are visited more?
m1 <- glmmTMB(prev_Stipes ~ Visit_interval + (1|Card_number),
              data = kelp)
summary(m1)

ggplot(kelp, aes(x = Visit_interval, y = prev_Stipes)) +
  geom_point()

#sst through seasons
m2 <- glmmTMB(sst_mean ~ Month + (1|Year) + (1|Card_number),
              data = kelp)
summary(m2)

ggplot(kelp, aes(x = Month, y = sst_mean)) +
  geom_abline(slope = 0.26, intercept = 14.86) +
  geom_point()

#wave_power through seasons
m3 <- glmmTMB(wave_p_mean ~ Month + (1|Year) + (1|Card_number),
              data = kelp)
summary(m3)

ggplot(kelp, aes(x = Month, y = wave_p_mean)) +
  geom_abline(slope = -0.32, intercept = 4.54) +
  geom_point()
