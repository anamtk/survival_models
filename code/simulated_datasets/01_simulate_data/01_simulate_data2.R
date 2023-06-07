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

source(here("code",
            "00_functions",
            "simulation_functions.R"))
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
  rownames_to_column(var = "ID") %>%
  pivot_longer(V1:V10,
               names_to = "interval",
               values_to = "x") %>%
  mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
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
  rownames_to_column(var = "ID") %>%
  pivot_longer(V1:V10,
               names_to = "interval",
               values_to = "x") %>%
  mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
  mutate(interval = as.numeric(interval))
# 2c: High variation ------------------------------------------------------

#create an empty matrix of the number of sites (rows)
# and number of intervals (columns)
x.int.high <- matrix(NA, nrow = n.indiv, ncol = n.int)

for(i in 1:n.indiv){
  for(j in 1:n.int){
    #this sd creates ~ 1 unit of variation
    x.int.high[i,j] <-  rnorm(1, mean = x.nom[i], sd = 1)
  } }

#make a long dataframe for merging with y data
x.int.high.df <- as.data.frame(x.int.high) %>%
  rownames_to_column(var = "ID") %>%
  pivot_longer(V1:V10,
               names_to = "interval",
               values_to = "x") %>%
  mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
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

ps_low <- (exp(b0 + b1*x.int.low))/(1 + exp(b0 + b1*x.int.low))

# b. Med var data ---------------------------------------------------------

ps_med<- (exp(b0 + b1*x.int.med))/(1 + exp(b0 + b1*x.int.med))

# c. high var -------------------------------------------------------------

ps_high<- (exp(b0 + b1*x.int.high))/(1 + exp(b0 + b1*x.int.high))


# Step 4: Simulate y data -------------------------------------------------

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
  left_join(x.int.high.df, by = c("ID", "interval")) %>%
  mutate(t = 1)

# Export datasets ---------------------------------------------------------

write.csv(y.low.all2, here("data_outputs",
                           "simulated",
                           "02_analysis_ready",
                           "low_var_interval_data.csv"))

write.csv(y.med.all2, here("data_outputs",
                           "simulated",
                           "02_analysis_ready",
                           "med_var_interval_data.csv"))

write.csv(y.high.all2, here("data_outputs",
                           "simulated",
                           "02_analysis_ready",
                           "high_var_interval_data.csv"))


# Full survey summarised data ---------------------------------------------

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

