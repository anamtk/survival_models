# Simulated datasets for survival mdoels
# Ana Miller-ter Kuile
# April 3, 2023

#This script simulates datasets to try in the different
#survival models


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")
## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
## And loading them
for(i in package.list){library(i, character.only = T)}
times



# Path forward - datasets needed ------------------------------------------

#1.	Many surveys (~5 each), unequal 1-0 (25%:75% end fates)(poor fit) – this is our dataset
#2.	Many surveys (~5 each), equal 1-0 (50:50 end result)(some poor fit)
#3.	Few surveys (~2 each), unequal 1-0 (25:75) (should fit okay)
#4.	Few surveys (~2 each), equal 1-0 (50:50) (should fit well)


# All dataset values ------------------------------------------------------

set.seed(1)

#number of times thing was surveyed
times <- round(rnorm(500, mean = 4))
tot <- sum(times) #the total intervals from that value
#continuous predictor covariate
predictor <- rpois(tot, lambda =10)


# Uneven 1-0 dataset simulation -------------------------------------------

#setting intercept >0 makes more 1s than 0s in dataset
xbu <- 3.1 + 0.2*scale(predictor, center = TRUE, scale = FALSE)
pu <- 1/(1 + exp(-xbu))
yu <- rbinom(n = tot, size = 1, prob = pu)
#how many 0's?
zerou <- length(yu[which(yu == 0)])
#107/500 ~ 21%
#107/2020 ~ 5%

# Even 1-0 dataset simulation ---------------------------------------------

#setting intercept to closer to 0 makes ~equal 1-0 dataset
xbe <- 2 + 0.2*scale(predictor, center = TRUE, scale = FALSE)
pe <- 1/(1 + exp(-xbe))
ye <- rbinom(n = tot, size = 1, prob = pe)
#how many 0's?
zeroe <- length(ye[which(ye == 0)])
#276/500 ~55%
#276/2020 ~14%

# Make into dataframes ----------------------------------------------------

id <- 1:500

meta <- as.data.frame(cbind(id = id,
                            times = times))

id2 <- as.data.frame(rep(id, times)) %>%
  rename("ID" = 'rep(id, times)') %>%
  group_by(ID) %>%
  mutate(interval = 1:n()) %>%
  ungroup()

id3 <- id2 %>%
  group_by(ID) %>%
  summarise(interval = max(interval))


# Uneven 0's dataframe ----------------------------------------------------

uneven <- as.data.frame(cbind(predictor = predictor,
                              y = yu)) %>%
  arrange(y) %>%
  dplyr::select(-y)

uneven2 <- as.data.frame(cbind(predictor = predictor,
                              y = yu)) %>%
  arrange(y) 

zero_IDu <- sample(id, size = zerou, replace = F)
zero_IDu

idu <- id3 %>%
  mutate(fate = case_when(ID %in% zero_IDu ~ 0,
                          TRUE ~ 1)) %>%
  full_join(id2, by = c("ID", "interval")) %>%
  replace_na(list(fate = 1)) %>%
  arrange(fate) %>%
  bind_cols(uneven)


# Even 0's dataframe ------------------------------------------------------


even <- as.data.frame(cbind(predictor = predictor,
                              y = ye)) %>%
  arrange(y) %>%
  dplyr::select(-y)

even2 <- as.data.frame(cbind(predictor = predictor,
                               y = ye)) %>%
  arrange(y) 

zero_IDe <- sample(id, size = zeroe, replace = F)
zero_IDe

ide <- id3 %>%
  mutate(fate = case_when(ID %in% zero_IDe ~ 0,
                          TRUE ~ 1)) %>%
  full_join(id2, by = c("ID", "interval")) %>%
  replace_na(list(fate = 1)) %>%
  arrange(fate) %>%
  bind_cols(even)


# Fewer intervals ---------------------------------------------------------

times2 <- ceiling(times/2)

id4 <- as.data.frame(rep(id, times2)) %>%
  rename("ID" = 'rep(id, times2)') %>%
  group_by(ID) %>%
  mutate(interval = 1:n()) %>%
  ungroup()

num <- length(id4$interval)
# Uneven dataset ----------------------------------------------------------


uneven3 <- uneven2[1:num, 1:2] %>%
  dplyr::select(-y)

idu2 <- id4 %>%
  mutate(fate = case_when(ID %in% zero_IDu ~ 0,
                          TRUE ~ 1)) %>%
  full_join(id4, by = c("ID", "interval")) %>%
  replace_na(list(fate = 1)) %>%
  arrange(fate) %>%
  bind_cols(uneven3)


# Even dataset ------------------------------------------------------------


even3 <- even2[1:num, 1:2] %>%
  dplyr::select(-y)

ide2 <- id4 %>%
  mutate(fate = case_when(ID %in% zero_IDe ~ 0,
                          TRUE ~ 1)) %>%
  full_join(id4, by = c("ID", "interval")) %>%
  replace_na(list(fate = 1)) %>%
  arrange(fate) %>%
  bind_cols(even3)


