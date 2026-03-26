
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "readxl", 'glmmTMB',
                  'lubridate')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------


climate <- read.csv(here("data_outputs",
                         "01_whwonests",
                         "01_cleaning",
                         "03_dailyWeather_nests.csv"))


climate %>%
ggplot(aes(x = as.Date(Date), y = tmax, group = Nest_ID)) +
  geom_line() +
  facet_wrap(~Year_located, scales = "free_x")

climate2 <- climate %>%
  mutate(julian = yday(as.Date(Date))) %>%
  mutate(tmax = scale(tmax))



m1 <- glmmTMB(tmax ~ julian + (1|Year_located),
              data = climate2)

summary(m1)

sum <- summary(m1)

coef <- as.data.frame(sum$coefficients$cond) %>%
  dplyr::select(Estimate)

b0 <- coef[1,]

re <- m1$fit$parfull
re_names <- names(re)

re_df <- bind_cols(re, re_names) %>%
  rename("value" = "...1",
         "re" = "...2") %>%
  filter(!str_detect(re, "beta|theta")) %>%
  mutate(Year_located = 2012:(2011+n())) %>%
  dplyr::select(-re) %>%
  rename("re" = "value") %>%
  mutate(b0 = coef[1,],
         b = coef[2,])

days <- climate2 %>%
  group_by(Nest_ID, Year_located) %>%
  reframe(min = min(julian),
          max = max(julian)) %>%
  ungroup() %>%
  group_by(Year_located) %>%
  reframe(min = min(min),
          max = max(max)) %>%
  ungroup()

all_df <- days %>%
  left_join(re_df, by = "Year_located") %>%
  rowwise() %>%
  mutate(min_val = b*min + b0 + re,
         max_val = b*max + b0 + re) %>%
  ungroup() %>%
  reframe(min_temp = mean(min_val),
          max_temp = mean(max_val))
