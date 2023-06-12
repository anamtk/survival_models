# Simulating datasets for testing models
# Ana Miller-ter Kuile
# May 17, 2023

#this script simulates datasets for testing the different
#survival model structures

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#I put this function into the code below
# source(here("code",
#             "00_functions",
#             "simulation_functions.R"))

# Data specifications -----------------------------------------------------

# -Response data need to be 1-0 data for a set of individuals
## across a certain number of survey sub-intervals
# -Response is dependent on a covariate 
## (could also have strong vs. weak covariate effects as
## an additional level of variation later)
# - we have a set of individuals that we'll be modeling for
## there is ~large variation *across* individuals in the 
## nominal covariate
# - We will vary the amount of variation *within individuals* 
## *across* visit intervals

#STEPS:
#DO JUST ONCE:
# 1. Simulate individual level "nominal" x:
## X[i] ~ N(0,1)
# 2. Then simulate interval level x given the site
## "nominal" value:
### a. small variation: X[i,j] ~ N(X[i], ?) #small, 0.01?
### b. med variation: X[i,j] ~ N(X[i], ?) #+-1 unit of X
### c. large variation (comparable to across site)
### X[i,j] ~ N(X[i], 1)
# 3. plug X[i,j]'s into p_s logistic regression equation
#REPEAT FOR 100 "DATASETS":
# 4. Generate y's via: y[i,j] ~ bern(p_s[i,j])
### a. create each individual data timeseries from generating
### these y's
### b. cut off all series after the first 0
### OR after the cutoff of the max number of intervals

# Step 1: Create nominal X ------------------------------------------------

#number of individuals to begin with
n.indiv <- 300
#number of intervals per individual
n.int <- 10

#simulate mean x per site for each individual
x.nom <- rnorm(n.indiv, mean = 0, sd = 1)

# Step 2: Simulate interval x's -------------------------------------------

# 2a: low variability -----------------------------------------------------

#create an empty matrix of the number of sites (rows)
# and number of intervals (columns)
x.int.low <- matrix(NA, nrow = n.indiv, ncol = n.int)

for(i in 1:n.indiv){
  for(j in 1:n.int){
    #this sd creates a 0.1 level variation 
    x.int.low[i,j] <-  rnorm(1, mean = x.nom[i], sd = 0.01)
  } }

#make this a long dataframe to merge later
x.int.low.df <- as.data.frame(x.int.low) %>%
  #make each rowname an ID
  rownames_to_column(var = "ID") %>%
  #make a long dataframe
  pivot_longer(V1:V10,
               names_to = "interval",
               values_to = "x") %>%
  #remove the "V" from every interval name
  mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
  #make interval a numeric value
  mutate(interval = as.numeric(interval))

# 2b: med variability -----------------------------------------------------

#create an empty matrix of the number of sites (rows)
# and number of intervals (columns)
x.int.med <- matrix(NA, nrow = n.indiv, ncol = n.int)

for(i in 1:n.indiv){
  for(j in 1:n.int){
    #this sd creates ~ 1 unit of variation
    x.int.med[i,j] <-  rnorm(1, mean = x.nom[i], sd = 0.5)
  } }

#make a long dataframe for merging with y data
x.int.med.df <- as.data.frame(x.int.med) %>%
  #make rownames an "ID" column
  rownames_to_column(var = "ID") %>%
  #make a long "tidy" dataframe 
  pivot_longer(V1:V10,
               names_to = "interval",
               values_to = "x") %>%
  #remove "V" in intervals
  mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
  #make intervals numeric
  mutate(interval = as.numeric(interval))

# 2c: High variation ------------------------------------------------------

#create an empty matrix of the number of sites (rows)
# and number of intervals (columns)
x.int.high <- matrix(NA, nrow = n.indiv, ncol = n.int)

for(i in 1:n.indiv){
  for(j in 1:n.int){
    #this sd creates variation similar to cross-site
    x.int.high[i,j] <-  rnorm(1, mean = x.nom[i], sd = 1)
  } }

#make a long dataframe for merging with y data
x.int.high.df <- as.data.frame(x.int.high) %>%
  #make an "ID" column from row number
  rownames_to_column(var = "ID") %>%
  #make dataframe long "tidy"
  pivot_longer(V1:V10,
               names_to = "interval",
               values_to = "x") %>%
  #get rid of "V" in inteval names
  mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
  #make interval numeric
  mutate(interval = as.numeric(interval))
# Step 3: plug x's into prob regression -----------------------------------

# Set betas ---------------------------------------------------------------

#we make betas for fairly high survival (usually
# above .5)

b0 <- 1.38 #gets mean value survival fairly high
b1 <- 0.5 #this value makes sure that ps is always >0.5 in the range of x

#down the road, could make this portion a stochastic
#process with betas centering around a value, and
#regenerating them 100 times as well, for now we'll
#make them deterministic though

# a. Low var data ---------------------------------------------------------

#Get the survival probability for each interval for each individual site
#for the low variability data
ps_low <- (exp(b0 + b1*x.int.low))/(1 + exp(b0 + b1*x.int.low))

# b. Med var data ---------------------------------------------------------

#Get the survival probability for each interval for each individual site
#for the med variability data
ps_med<- (exp(b0 + b1*x.int.med))/(1 + exp(b0 + b1*x.int.med))

# c. high var -------------------------------------------------------------

#Get the survival probability for each interval for each individual site
#for the high variability data
ps_high<- (exp(b0 + b1*x.int.high))/(1 + exp(b0 + b1*x.int.high))


# Step 4: Simulate y data -------------------------------------------------

#need to create a function that simulates a dataset of the 
#number of individuals and intervals per each with the probability
#matrix defined as one of the ones above (low, med, high)
y_function2 <- function(n.ind, 
                        n.t,
                        prob_matrix){
  
  #make an empty matrix for the ya data
  matrix <- matrix(NA, nrow = n.ind, ncol = n.t)
  
  #fill that matrix based on the probability matrix
  for(i in 1:n.ind){
    for(j in 1:n.t){ 
      matrix[i,j] <-  rbinom(1, size = 1, prob_matrix[i,j])
    } }
  
  #create a long format that we can filter
  df1 <- as.data.frame(matrix) %>%
    #make each row a numerical "ID"
    rownames_to_column(var = "ID") %>%
    #make a long "tidy" dataframe
    pivot_longer(V1:V10,
                 names_to = "interval",
                 values_to = "fate") %>%
    #remvoe the "V" from each interval
    mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
    #make interval numeric
    mutate(interval = as.numeric(interval))
  
  #DATA TRUNCATION STEP
  #the data are 1-0 data across a set of survey intervals per individual
  #Right now, this means that an indiviudal can be 1 - 0 - 1, which is not
  #possible (once dead, they stay dead)
  #these next steps get rid of any intervals after the first "dead" 
  #observation per individual
  
  #which individuals are alive at the end of the intervals?
  #subset just those individuals in a dataframe
  alive <- df1 %>%
    group_by(ID) %>%
    #filter out IDs where all intervals == 1
    filter(all(fate == 1)) %>%
    ungroup() 
  
  #which individuals are dead at some point in any
  #interval?
  #filter those out and then filter out any intervals
  #after the first one in which they are dead
  dead <- df1 %>%
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
  
  #combine these into one dataframe
  df2 <- alive %>%
    bind_rows(dead) 
  
  return(df2)
  
}

# a. low variation data ---------------------------------------------------

#use a custom function that creates 100 datatsets of y
#based on the probability of survival from the low variation dataset
#this function also removes any intervals after a 0 per individual
y.low <- lapply(1:100, #how many datasets
                function(i)
                {
                  y_function2(n.ind = n.indiv, #creates 1-0 for this many individauls,
                              n.t = n.int, #and this many intervals each
                              prob_matrix = ps_low) #based on this probability
                } )

#create this into one long df with "dataset" ID column corresponding
#to which element in the list above that dataframe is
y.low.all <- bind_rows(y.low, .id = "Dataset")

#bind that with the environmental data - such that there is only
#an X for a dataset for the intervals in which there are y data
y.low.all2 <- y.low.all %>%
  #combine y with x data
  left_join(x.int.low.df, by = c("ID", "interval")) %>%
  mutate(t = 1)

# Med variation data ------------------------------------------------------

#use a custom function that creates 100 datatsets of y
#based on the probability of survival from the med variation dataset
#this function also removes any intervals after a 0 per individual
y.med <- lapply(1:100, #how many datasets
                function(i)
                {
                  y_function2(n.ind = n.indiv, #creates 1-0 for this many individauls,
                              n.t = n.int, #and this many intervals each
                              prob_matrix = ps_med) #based on this probability
                } )

#create this into one long df with "dataset" ID column corresponding
#to which element in the list above that dataframe is
y.med.all <- bind_rows(y.med, .id = "Dataset")

#bind that with the environmental data - such that there is only
#an X for a dataset for the intervals in which there are y data
y.med.all2 <- y.med.all %>%
  #attach OG x data
  left_join(x.int.med.df, by = c("ID", "interval")) %>%
  mutate(t = 1)

# High variation data -----------------------------------------------------


#use a custom function that creates 100 datatsets of y
#based on the probability of survival from the high variation dataset
#this function also removes any intervals after a 0 per individual
y.high <- lapply(1:100, #how many datasets
                 function(i)
                 {
                   y_function2(n.ind = n.indiv, #creates 1-0 for this many individauls,
                               n.t = n.int, #and this many intervals each
                               prob_matrix = ps_high) #based on this probability
                 } )

#create this into one long df with "dataset" ID column corresponding
#to which element in the list above that dataframe is
y.high.all <- bind_rows(y.high, .id = "Dataset")

#bind that with the environmental data - such that there is only
#an X for a dataset for the intervals in which there are y data
y.high.all2 <- y.high.all %>%
  #attach OG x data
  left_join(x.int.high.df, by = c("ID", "interval")) %>%
  mutate(t = 1)

# Export datasets ---------------------------------------------------------

# write.csv(y.low.all2, here("data_outputs",
#                            "simulated",
#                            "02_analysis_ready",
#                            "low_var_interval_data.csv"))
# 
# write.csv(y.med.all2, here("data_outputs",
#                            "simulated",
#                            "02_analysis_ready",
#                            "med_var_interval_data.csv"))
# 
# write.csv(y.high.all2, here("data_outputs",
#                            "simulated",
#                            "02_analysis_ready",
#                            "high_var_interval_data.csv"))


# Full survey summarised data ---------------------------------------------

#for the total exposure model, need the above data to be summarised
#across the entire "lifetime" of the indiviudal, which is what
#these next steps do

low.var.tot <- y.low.all2 %>%
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
med.var.tot <- y.med.all2 %>%
  group_by(Dataset, ID) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         t= n(),
         x = mean(x)) %>%
  ungroup() %>%
  distinct(Dataset, ID, fate, t, x)

#do the same for the high variability dataset
high.var.tot <- y.high.all2 %>%
  group_by(Dataset, ID) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         t= n(),
         x = mean(x)) %>%
  ungroup() %>%
  distinct(Dataset, ID, fate, t, x)

# Export total exposure datasets ------------------------------------------

# write.csv(low.var.tot, here("data_outputs",
#                            "simulated",
#                            "02_analysis_ready",
#                            "low_var_total_data.csv"))
# 
# write.csv(med.var.tot, here("data_outputs",
#                            "simulated",
#                            "02_analysis_ready",
#                            "med_var_total_data.csv"))
# 
# write.csv(high.var.tot, here("data_outputs",
#                             "simulated",
#                             "02_analysis_ready",
#                             "high_var_total_data.csv"))

