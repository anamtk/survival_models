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
#of length of total number of intervals
predictor <- rpois(tot, lambda =10)


# Uneven 1-0 dataset simulation -------------------------------------------

#setting intercept >0 makes more 1s than 0s in dataset
#get the linear model
xbu <- 3.1 + 0.2*scale(predictor, center = TRUE, scale = FALSE)
#transform to logistic 
pu <- 1/(1 + exp(-xbu))
#get y values from that
yu <- rbinom(n = tot, size = 1, prob = pu)
#how many 0's?
zerou <- length(yu[which(yu == 0)])
#107/500 ~ 21%
#107/2020 ~ 5%

# Even 1-0 dataset simulation ---------------------------------------------

#setting intercept to closer to 0 makes ~equal 1-0 dataset
#get the linear model
xbe <- 2 + 0.2*scale(predictor, center = TRUE, scale = FALSE)
#transform to logistic
pe <- 1/(1 + exp(-xbe))
#get y's
ye <- rbinom(n = tot, size = 1, prob = pe)
#how many 0's?
zeroe <- length(ye[which(ye == 0)])
#276/500 ~55%
#276/2020 ~14%

# Make into dataframes ----------------------------------------------------

#ID is the individual being tracked
id <- 1:500

#get info on how many times each individual was surveyed
meta <- as.data.frame(cbind(id = id,
                            times = times))

#get an id dataframe that has the id replicated for the 
#number of survey intervals per individual,
#then, give these intervals values 1:n()
id2 <- as.data.frame(rep(id, times)) %>%
  rename("ID" = 'rep(id, times)') %>%
  group_by(ID) %>%
  mutate(interval = 1:n()) %>%
  ungroup()

#get the maximum interval for each individual being tracked
id3 <- id2 %>%
  group_by(ID) %>%
  summarise(interval = max(interval))


# Uneven 0's dataframe ----------------------------------------------------

#match the predictor and y values in a dataframe
uneven <- as.data.frame(cbind(predictor = predictor,
                              y = yu)) %>%
  #arrange so all 0 y's are first
  arrange(y) %>%
  #femove the y so we can bind it back up later
  dplyr::select(-y)

#this is just to double check that this is correct order for later
uneven2 <- as.data.frame(cbind(predictor = predictor,
                              y = yu)) %>%
  arrange(y) 

#get a random set of individuals whose final fate was a 0
zero_IDu <- sample(id, size = zerou, replace = F)
zero_IDu #which numbers are they

#combine for final dataset
dat_uneven <- id3 %>% #the max interval per individual dataset
  #when the ID is in the final fate = 0, set to 0, otherwise to 1
  mutate(fate = case_when(ID %in% zero_IDu ~ 0,
                          TRUE ~ 1)) %>%
  #combine with the full survey dataset, so that the final outcome is
  #tracked
  full_join(id2, by = c("ID", "interval")) %>%
  #set all previous outcomes to 1
  replace_na(list(fate = 1)) %>%
  #arrange by 0's to 1's
  arrange(fate) %>%
  #bind to the predictor dataset in the same order
  bind_cols(uneven)


# Even 0's dataframe ------------------------------------------------------

#Do the same for the even number of 1-0 final dataset

even <- as.data.frame(cbind(predictor = predictor,
                              y = ye)) %>%
  arrange(y) %>%
  dplyr::select(-y)

even2 <- as.data.frame(cbind(predictor = predictor,
                               y = ye)) %>%
  arrange(y) 

zero_IDe <- sample(id, size = zeroe, replace = F)
zero_IDe

dat_even <- id3 %>%
  mutate(fate = case_when(ID %in% zero_IDe ~ 0,
                          TRUE ~ 1)) %>%
  full_join(id2, by = c("ID", "interval")) %>%
  replace_na(list(fate = 1)) %>%
  arrange(fate) %>%
  bind_cols(even)


# Fewer intervals ---------------------------------------------------------

#reduce the number of intervals per individual
times2 <- ceiling(times/2)

#create a dataframe that gets a row for each intdividual for each interval
id4 <- as.data.frame(rep(id, times2)) %>%
  rename("ID" = 'rep(id, times2)') %>%
  group_by(ID) %>%
  #give those intervals values 1:n()
  mutate(interval = 1:n()) %>%
  ungroup()

#get the total length of this dataframe for later
num <- length(id4$interval)

#get the maximum interval for each individual being tracked
id5 <- id4 %>%
  group_by(ID) %>%
  summarise(interval = max(interval))

# Uneven dataset ----------------------------------------------------------

#select the first columns of the uneven dataset, which encapsulates
# all the zeros (since it's an arranged dataframe), and now has the 
#same number of rows as the reduced dataframe
uneven3 <- uneven2[1:num, 1:2] %>%
  dplyr::select(-y)

#now do the same thing of comining the final fate list with
#info about who died, then join with the rest of the survey
#interval dataset, set NA to 1 fates, and bind with the 
#predictor vector
dat_uneven2 <- id5 %>%
  mutate(fate = case_when(ID %in% zero_IDu ~ 0,
                          TRUE ~ 1)) %>%
  full_join(id4, by = c("ID", "interval")) %>%
  replace_na(list(fate = 1)) %>%
  arrange(fate) %>%
  bind_cols(uneven3)


# Even dataset ------------------------------------------------------------


even3 <- even2[1:num, 1:2] %>%
  dplyr::select(-y)

dat_even2 <- id5 %>%
  mutate(fate = case_when(ID %in% zero_IDe ~ 0,
                          TRUE ~ 1)) %>%
  full_join(id4, by = c("ID", "interval")) %>%
  replace_na(list(fate = 1)) %>%
  arrange(fate) %>%
  bind_cols(even3)


# Get full-survey summaries per individual --------------------------------


# Uneven ------------------------------------------------------------------

#uneven full
dat_u <- id3 %>% #the max interval per individual dataset
  #when the ID is in the final fate = 0, set to 0, otherwise to 1
  mutate(fate = case_when(ID %in% zero_IDu ~ 0,
                          TRUE ~ 1)) %>%
  dplyr::select(-interval)

dat_uneven3 <- dat_uneven %>% 
  ungroup() %>%
  dplyr::select(ID, predictor) %>%
  group_by(ID) %>%
  summarise(predictor = mean(predictor)) %>%
  ungroup() %>%
  left_join(dat_u, by = "ID")

#uneven reduced
dat_u2 <- id5 %>% #the max interval per individual dataset
  #when the ID is in the final fate = 0, set to 0, otherwise to 1
  mutate(fate = case_when(ID %in% zero_IDu ~ 0,
                          TRUE ~ 1)) %>%
  dplyr::select(-interval)
 
dat_uneven4 <- dat_uneven2 %>% 
  ungroup() %>%
  dplyr::select(ID, predictor) %>%
  group_by(ID) %>%
  summarise(predictor = mean(predictor)) %>%
  ungroup() %>%
  left_join(dat_u2, by = "ID") 


# Even --------------------------------------------------------------------

#even full
dat_e <- id3 %>% #the max interval per individual dataset
  #when the ID is in the final fate = 0, set to 0, otherwise to 1
  mutate(fate = case_when(ID %in% zero_IDe ~ 0,
                          TRUE ~ 1)) %>%
  dplyr::select(-interval)

dat_even3 <- dat_even %>% 
  ungroup() %>%
  dplyr::select(ID, predictor) %>%
  group_by(ID) %>%
  summarise(predictor = mean(predictor)) %>%
  ungroup() %>%
  left_join(dat_e, by = "ID")

#even reduced
dat_e2 <- id5 %>% #the max interval per individual dataset
  #when the ID is in the final fate = 0, set to 0, otherwise to 1
  mutate(fate = case_when(ID %in% zero_IDe ~ 0,
                          TRUE ~ 1)) %>%
  dplyr::select(-interval)

dat_even4 <- dat_even2 %>% 
  ungroup() %>%
  dplyr::select(ID, predictor) %>%
  group_by(ID) %>%
  summarise(predictor = mean(predictor)) %>%
  ungroup() %>%
  left_join(dat_e2, by = "ID") 

# Export the final datasets -----------------------------------------------

write.csv(dat_uneven,here("data_outputs",
               "02_analysis_ready",
               "simulated",
               "more_ints_uneven.csv"))


write.csv(dat_uneven2,here("data_outputs",
               "02_analysis_ready",
               "simulated",
               "few_ints_uneven.csv"))

write.csv(dat_uneven3,here("data_outputs",
                           "02_analysis_ready",
                           "simulated",
                           "more_ints_uneven_total.csv"))

write.csv(dat_uneven4,here("data_outputs",
                           "02_analysis_ready",
                           "simulated",
                           "few_ints_uneven_total.csv"))


write.csv(dat_even, here("data_outputs",
               "02_analysis_ready",
               "simulated",
               "more_ints_even.csv"))


write.csv(dat_even2, here("data_outputs",
               "02_analysis_ready",
               "simulated",
               "few_ints_even.csv"))

write.csv(dat_even3, here("data_outputs",
                          "02_analysis_ready",
                          "simulated",
                          "more_ints_even_total.csv"))

write.csv(dat_even4, here("data_outputs",
                          "02_analysis_ready",
                          "simulated",
                          "few_ints_even_total.csv"))
