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
            "tidy_functions.R"))

# Load data ---------------------------------------------------------------

low <- read.csv(here("data_outputs",
                     "simulated",
                     "02_analysis_ready",
                     "lowvar_data.csv"))

# Prep data objects for models --------------------------------------------


# Model 1 -----------------------------------------------------------------

n.indiv <- low %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

y <- low %>%
  distinct(ID, y_end) %>%
  dplyr::select(y_end) %>%
  as_vector()

x1 <- low %>%
  group_by(ID) %>%
  summarise(x1_sc = mean(x1_sc)) %>%
  dplyr::select(x1_sc) %>%
  as_vector()

t <- low %>%
  group_by(ID) %>%
  summarise(t = sum(t)) %>%
  dplyr::select(t) %>%
  as_vector()

data1 <- list(n.indiv = n.indiv,
             y = y,
             x1 = x1,
             t = t)

saveRDS(data1, here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model1low_JAGS_data.RDS"))


# Model 2 -----------------------------------------------------------------

n.indiv <- low %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()

n.t <- low %>%
  distinct(ID, n.t) %>%
  dplyr::select(n.t) %>%
  as_vector()
#matrix
y <- low %>%
  dplyr::select(ID, int, y) %>%
  pivot_wider(names_from = int,
              values_from = y) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

#matrix
t <- low %>%
  dplyr::select(ID, int, t) %>%
  pivot_wider(names_from = int,
              values_from = t) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

#matrix
x1 <- low %>%
  dplyr::select(ID, int, x1_sc) %>%
  pivot_wider(names_from = int,
              values_from = x1_sc) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

data2 <- list(n.indiv = n.indiv,
             n.t = n.t,
             y = y,
             x1 = x1,
             t = t)

saveRDS(data2, here("data_outputs",
                   "simulated",
                   "03_jags_input_data",
                   "Model2low_JAGS_data.RDS"))


# Model 3 -----------------------------------------------------------------

low1 <- low %>%
  arrange(n.t)

#numbers
n.indiv1 <- low1 %>%
  distinct(ID, n.t) %>%
  filter(n.t == 1) %>%
  tally() %>%
  as_vector()

n.indiv <- low1 %>%
  distinct(ID) %>%
  tally() %>%
  as_vector()
  
#vectors
n.t <- low1 %>%
  distinct(ID, n.t) %>%
  dplyr::select(n.t) %>%
  as_vector()

y <- low1 %>%
  distinct(ID, y_end) %>%
  column_to_rownames(var = "ID") %>%
  dplyr::select(y_end) %>%
  as_vector()

#matrices:
t <- low1 %>%
  dplyr::select(ID, int, t) %>%
  pivot_wider(names_from = int,
              values_from = t) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

x1 <- low1 %>%
  dplyr::select(ID, int, x1_sc) %>%
  pivot_wider(names_from = int,
              values_from = x1_sc) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

data3 <- list(n.indiv = n.indiv,
              n.indiv1 = n.indiv1,
             n.t = n.t,
             y = y,
             x1 = x1,
             t = t)

saveRDS(data3, here("data_outputs",
                    "simulated",
                    "03_jags_input_data",
                    "Model3low_JAGS_data.RDS"))

