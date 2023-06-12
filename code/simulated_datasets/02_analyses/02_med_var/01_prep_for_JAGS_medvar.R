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

#interval data
med <- read.csv(here("data_outputs",
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
n.indiv <- med %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

#matrix of number of intervals per individual per dataset
n.t <- med %>%
  group_by(ID, Dataset) %>%
  tally(name = "n.t") %>%
  pivot_wider(names_from = Dataset,
              values_from = n.t) %>%
  arrange(ID) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix

#Creating arrays for other data bits
n.ind <- 300 #individuals
n.int <- 10 #max number of intervals per individual
n.data <- 100 #number of datasets

ID <- med$ID #define a vector of ids
interval <- med$interval #vector of intervals
Dataset <- med$Dataset #vector of datasets

#create empty y array with the appropriate dimensions
y <- array(NA, dim = c(n.ind, n.int, n.data))

#fill array based on id, interval, and dataset for each row
for(i in 1:dim(med)[1]){
  y[ID[i], interval[i], Dataset[i]] <- med[i, 5] #column 5 is y data
}

#Do the same for the t data array
t <- array(NA, dim = c(n.ind, n.int, n.data))

for(i in 1:dim(med)[1]){
  t[ID[i], interval[i], Dataset[i]] <- med[i, 7]
}

#do the same as for the y data as for the x data array
x <- array(NA, dim = c(n.ind, n.int, n.data))

for(i in 1:dim(med)[1]){
  x[ID[i], interval[i], Dataset[i]] <- med[i, 6]
}

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

#number of individuals per dataset with only one interval
n.indiv1 <- med %>%
  group_by(ID, Dataset) %>%
  mutate(n.t = n()) %>%
  distinct(ID, Dataset, n.t) %>%
  filter(n.t == 1) %>%
  group_by(Dataset) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()

#total individuals in datasets
n.indiv <- med %>%
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

n.data <- 100 #number of datasets

#empty list for the number of surveys per individual
n.int <- list()

#fill with custom function that filters out each dataset at a 
#time and gets this vector, appending it into our list
for(i in 1:n.data){
  n.int[[i]] <- n.t_fun(df = med, dataset = i)
}

#bind this list of vectors together into a matrix
n.t <- bind_cols(n.int) %>%
  as.matrix()

#empty list for covariates for individuals x survey periods 
# x datasets array
xs <- list()

#use custom x_fun - which exports a matrix per dataset
for(i in 1:n.data){
  xs[[i]] <- x_fun(df = med, dataset = i)
}

#make this into an array with the correct dimensions
x <- array(unlist(xs), dim = c(nrow(xs[[1]]), ncol(xs[[1]]), length(xs)))

#repeat the x process to get the ID x interval x Dataset array for t
ts <- list()

#use t_fun for this one
for(i in 1:n.data){
  ts[[i]] <- t_fun(df = med, dataset = i)
}

t <- array(unlist(ts), dim = c(nrow(ts[[1]]), ncol(ts[[1]]), length(ts)))

#use y_fun to get a matrix of y x Datasets
ys <- list()

for(i in 1:n.data){
  ys[[i]] <- y_fun(df = med, dataset = i)
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
              t = t)

#export it
saveRDS(data3, here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model3med_JAGS_data.RDS"))

