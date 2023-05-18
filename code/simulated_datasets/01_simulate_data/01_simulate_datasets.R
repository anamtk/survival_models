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
            "tidy_functions.R"))


# Data specifications -----------------------------------------------------

# -Response data need to be 1-0 data for a set of individuals
## across a certain number of survey sub-intervals
# -Response is dependent on two covariates and their interaction
## (right now, I just have one covariate as a test)
# -These covariates vary at three different levels for the 
## three different datasets (low, med, high)

# Simulate data -----------------------------------------------------------

set.seed(1)

# Create structure of individuals-visits ----------------------------------

#number of individuals in dataset
n.indiv <- 300
#number of times each individual was smpaled, with max 4
n.t <- sample(4, 300, replace = T)

#get a metadata file of the total number of surveys/individual
meta <- as.data.frame(cbind(ID = 1:300,
                            n.t = n.t))


# Generate predictor variables --------------------------------------------

#predictor variable, t
x1 <- rnorm(300, mean = 23, sd = 6)

ts <- c(1:10)
#bind predictor with other data, scale predictor
meta <- meta %>%
  bind_cols(x1= x1) %>%
  mutate(x1_sc = scale(x1)) %>%
  rowwise() %>%
  #get number of days/time steps per interval
  mutate(t = sample(ts,size = 1, replace = T)) 


# Simulate final fate data ------------------------------------------------

#final fate data simulation from logistic regression
intercept <- 6.81
beta <- 4.85 #effect of scaled t

#get regression model defined
z <-  intercept + beta*(meta$x1_sc)       # linear combination with a bias
#inverse logit that model to get survival probability
pr <-  (exp(z)/(1+exp(z)))^meta$t        # pass through an inv-logit function

#bind the probabilities and then calcluate the
#end fates
meta2 <- meta %>%
  bind_cols(pr = pr) %>%
  mutate(y_end = rbinom(1, 1, pr))


#get full visit dataset by repeating id by the n.t for the ID
visits <- as.data.frame(rep(meta$ID, meta$n.t)) %>%
  rename("ID" = 'rep(meta$ID, meta$n.t)')%>%
  left_join(meta2, by = "ID") %>%
  group_by(ID) %>%
  #get a 1:n interval number per ID
  mutate(int = 1:n()) %>%
  #set y to y_end when last interval
  mutate(y = case_when(int == max(int) ~ y_end)) %>%
  ungroup() %>%
  #set all other intervals for y = 1
  replace_na(list(y = 1)) %>%
  group_by(ID) %>%
  #set NA all values we don't actually know for
  #response data 
  mutate(x1_sc = case_when(int == max(int) ~ x1_sc,
                          TRUE ~ NA_real_)) %>%
  ungroup() %>%
  dplyr::select(-pr, -x1) %>%
  rowwise() %>%
  mutate(x1_sc = as.numeric(x1_sc))

#what is the mean predictor by end fate category?
meta2 %>%
  group_by(y_end) %>%
  summarise(mean = mean(x1_sc, na.rm = T),
            min = min(x1_sc),
            max = max(x1_sc),
            sd = sd(x1_sc))


# Fill in predictor for all survey intervals ------------------------------

#get a dataframe with low variation in the predictor variable
low_vardf <- visits %>%
  rowwise() %>%
  #set t_sc for all y=1 intervals with NAs currently
  #to be around the mean of successes for end fates, with
  #sd that we can vary for three different datasets
  mutate(x1_sc = case_when(is.na(x1_sc) ~ rnorm(1, mean = 0.3, sd = 0.1),
                          TRUE ~ x1_sc))
  

#check this model
model_low <- glm(y ~ x1_sc,
             data = low_vardf,
             offset = low_vardf$t,
             family = "binomial")

summary(model_low)


med_vardf <- visits %>%
  rowwise() %>%
  #set t_sc for all y=1 intervals with NAs currently
  #to be around the mean of successes for end fates, with
  #sd that we can vary for three different datasets
  mutate(x1_sc = case_when(is.na(x1_sc) ~ rnorm(1, mean = 0.226, sd = .5),
                          TRUE ~ x1_sc))


model_med <- glm(y ~ x1_sc,
                 data = med_vardf,
                 offset = med_vardf$t,
                 family = "binomial")

summary(model_med)



high_vardf <- visits %>%
  rowwise() %>%
  #set t_sc for all y=1 intervals with NAs currently
  #to be around the mean of successes for end fates, with
  #sd that we can vary for three different datasets
  mutate(x1_sc = case_when(is.na(x1_sc) ~ rnorm(1, mean = 0.226, sd = 1),
                          TRUE ~ x1_sc))


model_high <- glm(y ~ x1_sc,
                 data = high_vardf,
                 offset = high_vardf$t,
                 family = "binomial")

summary(model_high)

# Export final datasets ---------------------------------------------------

write.csv(low_vardf, 
          here("data_outputs",
               "simulated",
               "lowvar_data.csv"))

write.csv(med_vardf,
          here("data_outputs",
               "simulated",
               "medvar_data.csv"))

write.csv(high_vardf,
          here("data_outputs",
               "simulated",
               "highvar_data.csv"))


#END SCRIPT
