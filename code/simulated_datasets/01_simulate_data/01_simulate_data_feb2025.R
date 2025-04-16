
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

set.seed(1)
# Data specifications -----------------------------------------------------

#1. generate time series of covariate for 3 sites from zero normal
#2. standardize those time series
#3. multiply those time series by 
## small: 0.01
## medium: 0.1
## don’t have to do anything in the large variability case
#4. simulate survival of 100 individuals at each site with 
## a regression with b1 and b0 (see notes below)
#5. Repeat for a number of "datasets"

#double check and tweak b0 and b1 to get:
# - ~1/2 of individuals should die
# - b0 = 0 should make ~1/2 individuals die
# - how many die in first versus later intervals? want a ~large
## subset to survive through first interval but then ~1/2
## to die by end
# - check that survival variability varies more with
## higher covariate variation

# Generate covariate data -------------------------------------------------

#create three time series of covariate data,
#standardize them, then create them to 
#have variable amounts of variation 

#number of intervals
n.int <- 45

#number of individuals per "site"
n.indiv <- 30

#set b0 and b1
b0 <- 1.75
b1 <- 3.5

sim_xy_function <- function(n.t,
                         var_level,
                         b0, 
                         b1,
                         n.ind){
  
  #get x for each interval
  x1 <- rnorm(n.t, mean = 0, sd = 1)
  #scale that 
  x2 <- scale(x1)
  #multiply by the level of variation we want
  #in our dataset (low 0.01, med 0.1, high 1)
  x.vec <- x2*var_level
  
  #make that a dataframe and get
  #survival probability based on regression
  #and b0 and b1 values defined
  x.df <- as.data.frame(x.vec) %>%
    rename(x = V1) %>%
    mutate(day = row_number()) %>%
    rowwise() %>%
    mutate(ps = exp(b0 + b1*x)/(1 + exp(b0 + b1*x))) %>%
    ungroup()
  
  #generate y data for each of those intervals
  #with binomial likelihood based on survival
  #probability
  ind.df <- expand.grid(1:n.ind, 1:n.t) %>%
    rename(ID = Var1,
           interval = Var2) %>%
    left_join(x.df, by = "interval") %>%
    rowwise() %>%
    mutate(fate = rbinom(n(), 1, ps))
  
  # #DATA TRUNCATION STEP
  # #the data are 1-0 data across a set of survey intervals per individual
  # #Right now, this means that an indiviudal can be 1 - 0 - 1, which is not
  # #possible (once dead, they stay dead)
  # #these next steps get rid of any intervals after the first "dead" 
  # #observation per individual
  # 
  # #which individuals are alive at the end of the intervals?
  # #subset just those individuals in a dataframe
  alive <- ind.df %>%
    group_by(ID) %>%
    #filter out IDs where all intervals == 1
    filter(all(fate == 1)) %>%
    ungroup()
  # 
  # #which individuals are dead at some point in any
  # #interval?
  # #filter those out and then filter out any intervals
  # #after the first one in which they are dead
  dead <- ind.df %>%
    group_by(ID) %>%
    #filter out just IDs where there is at least one zero
    filter(any(fate == 0)) %>%
    #Get a True-False column that finds all 0's per indiivudla site
    mutate(first_0 = fate == 0 & !duplicated(fate == 0)) %>%
    # Find the first row with `first_0` in each group, filter out rows after it
    filter(row_number() <= min(which(fate == 0 & first_0 == TRUE))) %>%
    ungroup() %>%
    #remove gropuing column
    dplyr::select(-first_0)
  # 
  # #combine these into one dataframe
  df2 <- alive %>%
    bind_rows(dead)
  # 
  return(df2)
  
}

#based on one datset, seems to work pretty well
#for getting survival proportions that we want (see
#above notes)
# test_low <- sim_xy_function(b0 = b0,
#                             b1 = b1,
#                             n.t = 10,
#                             var_level = 0.01,
#                             n.ind = 30)
# 
# #total fate ~ equal
# test_low %>%
#   group_by(ID) %>%
#   filter(interval == max(interval)) %>%
#   group_by(fate) %>%
#   tally()
# 
# #some live through first interval
# test_low %>%
#   filter(interval == 1) %>%
#   group_by(fate) %>%
#   tally()
# 
# test_low %>%
#   group_by(interval, fate) %>%
#   tally()
# 
# test_med <- sim_xy_function(b0 = b0,
#                             b1 = b1,
#                             n.t = 10,
#                             var_level = 0.1,
#                             n.ind = 30)
# 
# #most live
# test_med %>%
#   group_by(ID) %>%
#   filter(interval == max(interval)) %>%
#   group_by(fate) %>%
#   tally()
# 
# #some live through first interval
# test_med %>%
#   filter(interval == 1) %>%
#   group_by(fate) %>%
#   tally()
# 
# test_med %>%
#   group_by(interval, fate) %>%
#   tally()
# 
# test_high <- sim_xy_function(b0 = b0,
#                              b1 = b1,
#                              n.t = 10,
#                              var_level = 1,
#                              n.ind = 30)
# 
# #fewer live here
# test_high %>%
#   group_by(ID) %>%
#   filter(interval == max(interval)) %>%
#   group_by(fate) %>%
#   tally()
# 
# #~1/2 live through first interval
# test_high %>%
#   filter(interval == 1) %>%
#   group_by(fate) %>%
#   tally()
# 
# test_high %>%
#   group_by(interval, fate) %>%
#   tally()
# 

# Step 3 - generate a set of 100 datasets at each level ---------------------

#low variation
y.low <- lapply(1:100, #how many datasets
                function(i)
                {
                  sim_xy_function(b0 = b0,
                                  b1 = b1,
                                  n.t = n.int,
                                  var_level = 0.01,
                                  n.ind = n.indiv) #based on this probability
                } )

#create this into one long df with "dataset" ID column corresponding
#to which element in the list above that dataframe is
y.low.all <- bind_rows(y.low, .id = "Dataset") %>%
  mutate(t = 1)

# y.low.all %>%
#   group_by(ID, Dataset) %>%
#   filter(interval == max(interval)) %>%
#   group_by(Dataset, fate) %>%
#   tally() %>%
#   mutate(fate = as.factor(fate)) %>%
#   ggplot(aes(x = fate, y = n)) +
#   geom_boxplot()

#med variation
y.med <- lapply(1:100, #how many datasets
                function(i)
                {
                  sim_xy_function(b0 = b0,
                              b1 = b1,
                              n.t = n.int,
                              var_level = 0.1,
                              n.ind = n.indiv) #based on this probability
                } )

#create this into one long df with "dataset" ID column corresponding
#to which element in the list above that dataframe is
y.med.all <- bind_rows(y.med, .id = "Dataset") %>%
  mutate(t = 1)

# y.med.all %>%
#   group_by(ID, Dataset) %>%
#   filter(interval == max(interval)) %>%
#   group_by(Dataset, fate) %>%
#   tally() %>%
#   mutate(fate = as.factor(fate)) %>%
#   ggplot(aes(x = fate, y = n)) +
#   geom_boxplot()

#high
y.high <- lapply(1:100, #how many datasets
                function(i)
                {
                  sim_xy_function(b0 = b0,
                              b1 = b1,
                              n.t = n.int,
                              var_level = 1,
                              n.ind = n.indiv) #based on this probability
                } )

#create this into one long df with "dataset" ID column corresponding
#to which element in the list above that dataframe is
y.high.all <- bind_rows(y.high, .id = "Dataset") %>%
  mutate(t = 1)

# y.high.all %>%
#   group_by(ID, Dataset) %>%
#   filter(interval == max(interval)) %>%
#   group_by(Dataset, fate) %>%
#   tally() %>%
#   mutate(fate = as.factor(fate)) %>%
#   ggplot(aes(x = fate, y = n)) +
#   geom_boxplot()


# Explore datasets --------------------------------------------------------

c <- y.high.all %>%
  ggplot(aes(x = interval, y= x, group = Dataset)) +
  geom_line()

b <- y.med.all %>%
  ggplot(aes(x = interval, y= x, group = Dataset)) +
  geom_line()

a <- y.low.all %>%
  ggplot(aes(x = interval, y= x, group = Dataset)) +
  geom_line()

#a+ b+ c

d <- y.low.all %>%
  ggplot(aes(x = interval, y = ps, group = Dataset)) +
  geom_line()

e <- y.med.all %>%
  ggplot(aes(x = interval, y = ps, group = Dataset)) +
  geom_line()

f <- y.high.all %>%
  ggplot(aes(x = interval, y = ps, group = Dataset)) +
  geom_line()

(a+ b+ c)/
(d + e +f)


# interval <- tibble(interval = sort(sample(1:5, 
#                                           n.d, 
#                                           replace = T)))


# Export datasets ---------------------------------------------------------

write.csv(y.low.all, here("data_outputs",
                           "simulated",
                           "02_analysis_ready",
                           "low_var_interval_data.csv"))

write.csv(y.med.all, here("data_outputs",
                           "simulated",
                           "02_analysis_ready",
                           "med_var_interval_data.csv"))

write.csv(y.high.all, here("data_outputs",
                            "simulated",
                            "02_analysis_ready",
                            "high_var_interval_data.csv"))


# Full survey summarised data ---------------------------------------------

#for the total exposure model, need the above data to be summarised
#across the entire "lifetime" of the indiviudal, which is what
#these next steps do

low.var.tot <- y.low.all %>%
  #for each dataset and id in each dataset
  group_by(Dataset, ID) %>%
  #get fate to be 0 if any are 0, otherwise 1
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         #total number of "time" (equal to number of intervals)
         t= n(),
         #x is mean of the value of x for that individual
         x = mean(x)) %>%
  ungroup() %>%
  #select dataset, id, fate, t, x
  distinct(Dataset, ID, fate, t, x)

#do the same process for the medium variability dataset
med.var.tot <- y.med.all %>%
  group_by(Dataset, ID) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         t= n(),
         x = mean(x)) %>%
  ungroup() %>%
  distinct(Dataset, ID, fate, t, x)

#do the same for the high variability dataset
high.var.tot <- y.high.all %>%
  group_by(Dataset, ID) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         t= n(),
         x = mean(x)) %>%
  ungroup() %>%
  distinct(Dataset, ID, fate, t, x)

# Export total exposure datasets ------------------------------------------

write.csv(low.var.tot, here("data_outputs",
                            "simulated",
                            "02_analysis_ready",
                            "low_var_total_data.csv"))

write.csv(med.var.tot, here("data_outputs",
                            "simulated",
                            "02_analysis_ready",
                            "med_var_total_data.csv"))

write.csv(high.var.tot, here("data_outputs",
                             "simulated",
                             "02_analysis_ready",
                             "high_var_total_data.csv"))

