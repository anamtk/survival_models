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
med.day <- read.csv(here("data_outputs",
                         "simulated",
                         "02_analysis_ready",
                         "med_var_daily_data.csv"))
#interval data
med.int <- read.csv(here("data_outputs",
                     "simulated",
                     "02_analysis_ready",
                     "med_var_interval_data.csv"))

#total data summarised
med.tot <- read.csv(here("data_outputs",
                         "simulated",
                         "02_analysis_ready",
                         "med_var_total_data.csv"))

# Prep data objects for models --------------------------------------------


# Model 1 -----------------------------------------------------------------

n.indiv <- med.tot %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

#matrix of individuals by dataset
y <- med.tot %>%
  dplyr::select(Dataset, ID, fate) %>%
  pivot_wider(names_from = Dataset,
              values_from = fate) %>%
  arrange(ID) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

#matrix of individuals x dataset
x <- med.tot %>%
  dplyr::select(Dataset, ID, x) %>%
  pivot_wider(names_from = Dataset,
              values_from = x) %>%
  arrange(ID) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

#matrix of individuals by dataset
t <- med.tot %>%
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
                    "Model1med_JAGS_data.RDS"))


# Model 2 -----------------------------------------------------------------

#number individuals in dataset
n.indiv <- med.int %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

#matrix of number of intervals per individual per dataset
n.t <- med.int %>%
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
n.int <- max(med.int$interval) #max number of intervals per individual
n.data <- 100 #number of datasets

ID <- med.int$ID #define a vector of ids
interval <- med.int$interval #vector of intervals
Dataset <- med.int$Dataset #vector of datasets

#create empty y array with the appropriate dimensions
y <- array(NA, dim = c(n.ind, n.int, n.data))

#fill array based on id, interval, and dataset for each row
for(i in 1:dim(med.int)[1]){
  y[ID[i], interval[i], Dataset[i]] <- med.int[i, 5] #column 5 is y data
}

#Do the same for the t data array
t <- array(NA, dim = c(n.ind, n.int, n.data))

for(i in 1:dim(med.int)[1]){
  t[ID[i], interval[i], Dataset[i]] <- med.int[i, 6]
}

#x is by site and dataset
x <- med.int %>%
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
                    "Model2med_JAGS_data.RDS"))


# Model 3 -----------------------------------------------------------------

datasets <- as.data.frame(1:100) %>%
  rename('Dataset' = `1:100`) 

#get ids so individuals with one interval are first
meddayid <- med.day %>%
  group_by(ID, Dataset) %>%
  mutate(n.t = length(unique(interval))) %>%
  arrange(n.t) %>%
  distinct(Dataset, ID) %>%
  ungroup() %>%
  group_by(Dataset) %>%
  mutate(ID2 = 1:n())

med.day2 <- med.day  %>%
  group_by(ID, Dataset) %>%
  mutate(n.t = length(unique(interval))) %>%
  arrange(n.t) %>%
  ungroup() %>%
  left_join(meddayid, by = c("Dataset", "ID")) %>%
  dplyr::select(-ID) %>%
  rename(ID = ID2) %>%
  group_by(ID, Dataset, interval) %>%
  arrange(ID, Dataset, interval, day) %>%
  mutate(day.of.int = 1:n()) %>%
  ungroup() %>%
  as.data.frame() 

#number of individuals per dataset with only one interval
n.indiv1 <- med.day2 %>%
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
n.indiv <- med.day2 %>%
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
n.t <- med.day2 %>%
  distinct(Dataset, ID, interval) %>%
  group_by(Dataset, ID) %>%
  filter(interval == max(interval)) %>%
  ungroup() %>%
  pivot_wider(names_from = "Dataset",
              values_from = "interval") %>%
  dplyr::select(-ID) %>%
  as.matrix()

n.interval <- med.day2 %>%
  ungroup() %>%
  distinct(Dataset, interval) %>%
  group_by(Dataset) %>%
  filter(interval == max(interval)) %>%
  arrange(Dataset) %>%
  ungroup() %>%
  dplyr::select(interval) %>%
  as_vector()

#FOR ALL ARRAYS: Indexing dataframes AND 
#variables to loop through for for loops
med.day2 <- as.data.frame(med.day2)
n.ind <- 30 #individuals
n.int <- max(med.day2$interval) #max number of intervals per individual
n.dy <- max(med.day2$day)
n.data <- 100 #number of datasets

#arrays:
#i - individual
#j - interval
#d - dataset
#k - day
day.index.df <- med.day2 %>%
  group_by(ID, interval, Dataset) %>%
  mutate(start = min(day),
         end = max(day)) %>%
  ungroup() %>% 
  distinct(ID, interval, Dataset, start, end) %>%
  as.data.frame()

day.index.df2 <- day.index.df %>%
  distinct(interval, Dataset, start, end)

#Id <- day.index.df$ID #define a vector of ids
Int <- day.index.df$interval #vector of days in interval
Datset <- day.index.df$Dataset

#start.day[i,j,d] values from 1-24 (more like 3-22)
start.day <- matrix(NA, nrow = n.int, ncol = n.data)

#fill array based on id, interval, day, and dataset for each row
for(i in 1:dim(day.index.df2)[1]){
  start.day[Int[i], Datset[i]] <- day.index.df2[i, 3] 
}

#end.day[i,j,d] values from 1-24 (more like 3-24)
end.day <- matrix(NA, nrow = n.int, ncol = n.data)

#fill array based on id, interval, day, and dataset for each row
for(i in 1:dim(day.index.df2)[1]){
  end.day[Int[i], Datset[i]] <- day.index.df2[i, 4] 
}

#n.days[i,d]
n.days <- med.day2 %>%
  group_by(Dataset) %>%
  filter(day == max(day)) %>%
  distinct(Dataset, day) %>%
  arrange(Dataset) %>%
  ungroup() %>%
  dplyr::select(day) %>%
  as_vector()

#x[i,k,d]

#ID <- med.day2$ID #define a vector of ids
day <- med.day2$day #vector of days in interval
datset <- med.day2$Dataset

#create empty y array with the appropriate dimensions
x <- matrix(NA, nrow = n.dy, ncol = n.data)

#fill array based on id, interval, day, and dataset for each row
for(i in 1:dim(med.day2)[1]){
  x[day[i], datset[i]] <- med.day2[i, 4] 
}

#use y_fun to get a matrix of y x Datasets
ys <- list()

for(i in 1:n.data){
  ys[[i]] <- y_fun(df = med.day2, dataset = i)
}

y <- bind_cols(ys) %>%
  as.matrix()

#create a list of all the data
data3 <- list(n.datasets = 100,
              n.indiv = n.indiv,
              n.indiv1 = n.indiv1,
              n.t = n.t,
              n.interval = n.interval,
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
                    "Model3med_JAGS_data.RDS"))

