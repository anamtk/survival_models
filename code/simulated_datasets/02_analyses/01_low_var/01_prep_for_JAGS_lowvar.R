# Prep simulated data for JAGS
# Ana Miller-ter Kuile
# May 17, 2023

# prep data from simulated datasets for jags


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  "jagsUI",
                  "coda",
                  "mcmcplots")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("code",
            "00_functions",
            "simulation_functions.R"))

# Load data ---------------------------------------------------------------

#daily data
low.day <- read.csv(here("data_outputs",
                              "simulated",
                              "02_analysis_ready",
                              "low_var_daily_data.csv"))

#interval data
low.int <- read.csv(here("data_outputs",
                     "simulated",
                     "02_analysis_ready",
                     "low_var_interval_data.csv"))

#total data summarised
low.tot <- read.csv(here("data_outputs",
                         "simulated",
                         "02_analysis_ready",
                         "low_var_total_data.csv"))

# Prep data objects for models --------------------------------------------

# Model 1 -----------------------------------------------------------------

n.indiv <- low.tot %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

#matrix of individuals by dataset
y <- low.tot %>%
  dplyr::select(Dataset, ID, fate) %>%
  pivot_wider(names_from = Dataset,
              values_from = fate) %>%
  arrange(ID) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

#matrix of individuals x dataset
x <- low.tot %>%
  dplyr::select(Dataset, ID, x) %>%
  pivot_wider(names_from = Dataset,
              values_from = x) %>%
  arrange(ID) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

#matrix of individuals by dataset
t <- low.tot %>%
  dplyr::select(Dataset, ID, t) %>%
  pivot_wider(names_from = Dataset,
              values_from = t) %>%
  arrange(ID) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

#list these data for model
data1 <- list(n.datasets = 100,
              n.indiv = n.indiv,
             y = y,
             x = x,
             t = t)

saveRDS(data1, here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model1low_JAGS_data.RDS"))


# Model 2 -----------------------------------------------------------------
gc()
rm(list = 'n.indiv', 'y', 'x', 't')
#number individuals in dataset
n.indiv <- low.int %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

#matrix of number of intervals per individual per dataset
n.t <- low.int %>%
  distinct(ID, Dataset,interval) %>%
  group_by(ID, Dataset) %>%
  tally(name = "n.t") %>%
  pivot_wider(names_from = Dataset,
              values_from = n.t) %>%
  arrange(ID) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#Creating arrays for other data bits
n.ind <- 30 #individuals
n.int <- max(low.int$interval) #max number of intervals per individual
n.data <- 100 #number of datasets

ID <- low.int$ID #define a vector of ids
interval <- low.int$interval #vector of intervals
Dataset <- low.int$Dataset #vector of datasets

#create empty y array with the appropriate dimensions
y <- array(NA, dim = c(n.ind, n.int, n.data))

#fill array based on id, interval, and dataset for each row
for(i in 1:dim(low.int)[1]){
  y[ID[i], interval[i], Dataset[i]] <- low.int[i, 5] #column 5 is y data
}

#Do the same for the t data array
t <- array(NA, dim = c(n.ind, n.int, n.data))

for(i in 1:dim(low.int)[1]){
  t[ID[i], interval[i], Dataset[i]] <- low.int[i, 6]
}

#x is by site and dataset
x <- low.int %>%
  distinct(Dataset, interval, x) %>%
  pivot_wider(names_from = "Dataset",
              values_from = 'x') %>%
  column_to_rownames(var = "interval") %>%
  as.matrix()

#list of data for model 2
data2 <- list(n.datasets = 100,
              n.indiv = n.indiv,
             n.t = n.t,
             y = y,
             x = x,
             t = t)

#save it
saveRDS(data2, here("data_outputs",
                   "simulated",
                   "03_jags_input_data",
                   "Model2low_JAGS_data.RDS"))


# Model 3 -----------------------------------------------------------------
gc()
rm(list = 'n.indiv', 'y', 'x', 't')

datasets <- as.data.frame(1:100) %>%
  rename('Dataset' = `1:100`) 

#need to make sure this dataset is ordered by individuals with
#one interval first and multiple intervals next
lowdayid <- low.day %>%
  group_by(ID, Dataset) %>%
  mutate(n.t = length(unique(interval))) %>%
  arrange(n.t) %>%
  distinct(Dataset, ID) %>%
  ungroup() %>%
  group_by(Dataset) %>%
  mutate(ID2 = 1:n())

low.day2 <- low.day %>%
  group_by(ID, Dataset) %>%
  mutate(n.t = length(unique(interval))) %>%
  arrange(n.t) %>%
  ungroup() %>%
  left_join(lowdayid, by = c("Dataset", "ID")) %>%
  dplyr::select(-ID) %>%
  rename(ID = ID2) %>%
  group_by(ID, Dataset, interval) %>%
  arrange(ID, Dataset, interval, day) %>%
  mutate(day.of.int = 1:n()) %>%
  ungroup() %>%
  as.data.frame() 

#number of individuals per dataset with only one interval
n.indiv1 <- low.day2 %>%
  distinct(ID, Dataset, n.t) %>%
  filter(n.t == 1) %>%
  group_by(Dataset) %>%
  tally() %>%
  ungroup() %>%
  full_join(datasets, by = "Dataset") %>%
  replace_na(list(n = 0)) %>%
  arrange(Dataset) %>%
  dplyr::select(n) %>%
  as_vector()

#total individuals in datasets
n.indiv <- low.day2 %>%
  ungroup() %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

#to make matrices and arrays, need to use some custom functions
#so that the fact that different individuals in different datasets
#have one survey or more than one survey (e.g. individual 5
# could be surveyed once in datasets 5, 40, and 87 and more than
# once in other datasets - with the loops in the model, this will
#only work in the model if we make sure that each dataset is arranaged
# with the 1s first and then +1s (so ignoring actual individual IDs
#but maintaining the covariate info that generated the data))
#matrix of individual x dataset

#n.t[i,d]
n.t <- low.day2 %>%
  distinct(Dataset, ID, interval) %>%
  group_by(Dataset, ID) %>%
  filter(interval == max(interval)) %>%
  ungroup() %>%
  pivot_wider(names_from = "Dataset",
              values_from = "interval") %>%
  dplyr::select(-ID) %>%
  as.matrix()

#FOR ALL ARRAYS: Indexing dataframes AND 
#variables to loop through for for loops
low.day2 <- as.data.frame(low.day2)
n.ind <- 30 #individuals
n.int <- max(low.day2$interval) #max number of intervals per individual
n.dy <- max(low.day2$day)
n.data <- 100 #number of datasets

#arrays:
#i - individual
#j - interval
#d - dataset
#k - day
day.index.df <- low.day2 %>%
  group_by(ID, interval, Dataset) %>%
  mutate(start = min(day),
         end = max(day)) %>%
  ungroup() %>% 
  distinct(ID, interval, Dataset, start, end) %>%
  as.data.frame()

Id <- day.index.df$ID #define a vector of ids
Int <- day.index.df$interval #vector of days in interval
Datset <- day.index.df$Dataset

#start.day[i,j,d] values from 1-24 (more like 3-22)
start.day <- array(NA, dim = c(n.ind, n.int, n.data))

#fill array based on id, interval, day, and dataset for each row
for(i in 1:dim(day.index.df)[1]){
  start.day[Id[i], Int[i], Datset[i]] <- day.index.df[i, 4] 
}

#end.day[i,j,d] values from 1-24 (more like 3-24)
end.day <- array(NA, dim = c(n.ind, n.int, n.data))

#fill array based on id, interval, day, and dataset for each row
for(i in 1:dim(day.index.df)[1]){
  end.day[Id[i], Int[i], Datset[i]] <- day.index.df[i, 5] 
}

#n.days[i,d]
n.days <- low.day2 %>%
  group_by(ID, Dataset) %>%
  filter(day == max(day)) %>%
  dplyr::select(ID, Dataset, day) %>%
  ungroup() %>%
  pivot_wider(names_from = "Dataset",
              values_from = "day") %>%
  column_to_rownames(var= "ID") %>%
  as.matrix()

#n.days[i,j,d]
# days <- list()
# 
# for(i in 1:n.data){
#   days[[i]] <- n_day_fun(df = low.day2, dataset = i)
# }

# n.days <- array(unlist(days),
#                  dim = c(nrow(days[[1]]),
#                          ncol(days[[1]]),
#                          length(days)))

#x[i,k,d]

ID <- low.day2$ID #define a vector of ids
day <- low.day2$day #vector of days in interval
datset <- low.day2$Dataset

#create empty y array with the appropriate dimensions
x <- array(NA, dim = c(n.ind, n.dy, n.data))

#fill array based on id, interval, day, and dataset for each row
for(i in 1:dim(low.day2)[1]){
  x[ID[i], day[i], datset[i]] <- low.day2[i, 4] 
}

# #x[i,j,k,d]
# low.day2 <- as.data.frame(low.day2)
# n.ind <- 30 #individuals
# n.int <- max(low.day2$interval) #max number of intervals per individual
# n.dy <- max(low.day2$day.of.int)
# n.data <- 100 #number of datasets
# 
# ID <- low.day2$ID #define a vector of ids
# inter<- low.day2$interval #vector of intervals
# dayofint <- low.day2$day.of.int #vector of days in interval
# datset <- low.day2$Dataset
# 
# #create empty y array with the appropriate dimensions
# x <- array(NA, dim = c(n.ind, n.int, n.dy, n.data))
# 
# #fill array based on id, interval, day, and dataset for each row
# for(i in 1:dim(low.day2)[1]){
#   x[ID[i], inter[i], dayofint[i], datset[i]] <- low.day2[i, 4] 
# }
#   

#use y_fun to get a matrix of y x Datasets
ys <- list()

for(i in 1:n.data){
  ys[[i]] <- y_fun(df = low.day2, dataset = i)
}

y <- bind_cols(ys) %>%
  as.matrix()

#create a list of all the data
data3 <- list(n.datasets = 100,
              n.indiv = n.indiv,
              n.indiv1 = n.indiv1,
             n.t = n.t,
             y = y,
             x = x,
             n.days = n.days,
             start.day = start.day,
             end.day = end.day,
             y2 = y)

#export it
saveRDS(data3, here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model3low_JAGS_data.RDS"))

