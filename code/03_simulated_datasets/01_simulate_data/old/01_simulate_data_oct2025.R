
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

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

#Define two different interval "lengths" for each
#dataset
#short: 5 days/years
#long: 10 days/years

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

#number of intervals (now days)
n.int <- 45

#number of individuals per "site"
n.indiv <- 30

#set b0 and b1
b0 <- 2.8
#b0 <- 2.1
b1 <- 3.5

#number of simulated datasets
n.datasets <- 100

#simulated x data for every day of every interval
#in every dataset and KEEP them all

#Then generated y data for those x's - then, 
#"fill" 0's for intervals in which the individual "died"
#since this is the way the surveys would have been performed
#so for the last interval, we get all days in that interval, 
#and all are set to 0 - since we don't know which day in
#the total survey interval an individual died. 

#then combine:
#1. daily x data for each individual in each dataset
#in each day with interval IDs recognized
#2. interval level survival data with dataset, ID, and 
#interval ID info

sim_x_function <- function(n.t,
                           var_level,
                           b0, 
                           b1,
                           short.length,
                           long.length){
  
  #based on climate data from WHWO dataset,
  #creatng the x so that it gradually changes 
  #throughout the series
  days <- as.data.frame(155:(155+n.t))
  
  names(days) <- "day"
  
  days2 <- days %>%
    rowwise() %>%
    #this regression is from the tmax data from WHWO
    #with a regression of temp by julian day
    #also set julian days based on time of year
    #those nests were surveyed
    mutate(x1 = rnorm(n(), mean(0.0749171*day - 6.2250967), sd =1))
  
  #get x for each interval (day)
  #x1 <- rnorm(n.t, mean = 0, sd = 1)
  #scale that 
  x2 <- scale(days2$x1)
  #multiply by the level of variation we want
  #in our dataset (low 0.01, med 0.1, high 1)
  x.vec <- x2*var_level
  
  
  interval.long = sort(rep(1:(ceiling(n.t/long.length)), 
                           long.length))
  
  interval.long <- interval.long[1:n.t]
  #divide it into intervals of ~3.5 days each
  interval <- tibble(interval.short = sort(rep(1:(ceiling(n.t/short.length)), 
                                         short.length))) %>%
    mutate(day = 1:n()) %>%
    bind_cols(tibble(interval.long))
  
  #make that a dataframe and get
  #survival probability based on regression
  #and b0 and b1 values defined
  x.df <- as.data.frame(x.vec) %>%
    rename(x = V1) %>%
    mutate(day = row_number()) %>%
    rowwise() %>%
    mutate(ps = exp(b0 + b1*x)/(1 + exp(b0 + b1*x))) %>%
    ungroup() %>%
    left_join(interval, by = "day")
  
  return(x.df)
  
}


sim_y_function <- function(n.ind,
                           n.t,
                           n.data,
                           x.df){
  
  ind.df <- expand.grid(1:n.ind, 1:n.t, 1:n.data)%>%
    rename(ID = Var1,
           day = Var2,
           Dataset = Var3) %>%
    mutate(Dataset = as.character(Dataset)) %>%
    left_join(x.df, by = c("day", "Dataset")) %>%
    rowwise() %>%
    mutate(fate = rbinom(n(), 1, ps)) %>%
    ungroup()
  
  alive <- ind.df %>%
    group_by(ID, Dataset) %>%
    #filter out IDs where all intervals == 1
    filter(all(fate == 1)) %>%
    ungroup()
  
  dead <- ind.df %>%
    #arrange(day) %>%
    group_by(ID, Dataset) %>%
    #filter out just IDs where there is at least one zero
    filter(any(fate == 0)) %>%
    #Get a True-False column that finds all 0's per indiivudla site
    mutate(first_0 = fate == 0 & !duplicated(fate == 0)) %>%
    #instead of following code, to keep covariate data for 
    #interval summary model, just set fates after the first 0 to be equal to NA
    # mutate(fate = case_when(row_number() <= min(which(fate == 0 & first_0 == TRUE)) ~ fate,
    #                         TRUE ~ NA_real_)) %>%
    # Find the first row with `first_0` in each group, filter out rows after it
    filter(row_number() <= min(which(fate == 0 & first_0 == TRUE))) %>%
    ungroup() %>%
    #remove gropuing column
    dplyr::select(-first_0)
  
  df2 <- alive %>%
    bind_rows(dead) %>%
    dplyr::select(ID, day, Dataset, interval.long, 
                  interval.short, fate)
  
  ind.df2 <- ind.df %>%
    dplyr::select(-fate) %>%
    left_join(df2, by = c("ID", "day", "Dataset",
                          "interval.long", 'interval.short'))%>%
    group_by(ID, interval.short, Dataset) %>%
    mutate(fate.short = case_when(any(fate == 0) ~ 0,
                            TRUE ~ fate)) %>%
    ungroup() %>%
    group_by(ID, interval.long, Dataset) %>%
    mutate(fate.long = case_when(any(fate == 0) ~ 0,
                                  TRUE ~ fate)) %>%
    filter(!is.na(fate)) %>%
    ungroup() 
  
  return(ind.df2)
  
}

#Step 3 - generate a set of 100 datasets at each level ---------------------

set.seed(1)

#low variation
x.low <- lapply(1:100, #how many datasets
                function(i)
                {
                  sim_x_function(b0 = b0,
                                 b1 = b1,
                                 n.t = n.int,
                                 var_level = 0.01,
                                 short.length = 5,
                                 long.length = 10) #based on this probability
                } )


x.low.df <- bind_rows(x.low, .id = "Dataset") 
set.seed(1)
y.low <- sim_y_function(n.ind = n.indiv,
                        n.t = n.int,
                        n.data = n.datasets,
                        x.df = x.low.df)

set.seed(1)
#med variation
x.med <- lapply(1:100, #how many datasets
                function(i)
                {
                  sim_x_function(b0 = b0,
                                 b1 = b1,
                                 n.t = n.int,
                                 var_level = 0.1,
                                 short.length = 5,
                                 long.length = 10) #based on this probability
                } )


x.med.df <- bind_rows(x.med, .id = "Dataset") 

set.seed(1)
y.med <- sim_y_function(n.ind = n.indiv,
                        n.t = n.int,
                        n.data = n.datasets,
                        x.df = x.med.df)

set.seed(1)

#high
x.high <- lapply(1:100, #how many datasets
                function(i)
                {
                  sim_x_function(b0 = b0,
                                 b1 = b1,
                                 n.t = n.int,
                                 var_level = 0.5,
                                 short.length = 5,
                                 long.length = 10) #based on this probability
                } )


x.high.df <- bind_rows(x.high, .id = "Dataset") 

set.seed(1)
y.high <- sim_y_function(n.ind = n.indiv,
                        n.t = n.int,
                        n.data = n.datasets,
                        x.df = x.high.df)

# Explore datasets --------------------------------------------------------
theme_set(theme_bw())
y.high %>%
  ggplot(aes(x = day, y= x, group = Dataset)) +
  geom_line()

y.med %>%
  ggplot(aes(x = day, y= x, group = Dataset)) +
  geom_line()

y.low %>%
  ggplot(aes(x = day, y= x, group = Dataset)) +
  geom_line()

# Truncate and add intervals ----------------------------------------------

#for the high var, the max survival is to:
(max_day <- y.high %>%
  filter(!is.na(fate)) %>%
  filter(day == max(day)) %>%
   distinct(day) %>%
   dplyr::select(day) %>%
   as_vector())

#set manually to try to get more "live guys"
max_day <- 25

(max_interval_short <- y.high %>%
  filter(!is.na(fate)) %>%
  filter(day == max(day)) %>%
  distinct(interval.short) %>%
  dplyr::select(interval.short) %>%
  as_vector() )

(max_interval_long <- y.high %>%
    filter(!is.na(fate)) %>%
    filter(day == max(day)) %>%
    distinct(interval.long) %>%
    dplyr::select(interval.long) %>%
    as_vector() )

daily_survival_function <- function(df,
                                    interval_name,
                                    max_interval_val,
                                    other_int_name,
                                    other_int_fate,
                                    interval_fate){
  
  df2 <- df %>%
    filter({{interval_name}} <= max_interval_val) %>%
    filter(day <= max_day) %>%
    dplyr::select(-fate,
                  -{{other_int_name}},
                  -{{other_int_fate}}
                  ) %>%
    rename('interval' = {{interval_name}},
           'fate' = {{interval_fate}})%>%
    filter(!is.na(fate))
  
  return(df2)

}

y.low.daily.short <- daily_survival_function(df = y.low,
                        interval_name = interval.short,
                        max_interval_val = max_interval_short,
                        other_int_name = 'interval.long',
                        other_int_fate = 'fate.long',
                        interval_fate = 'fate.short')

y.low.daily.long <- daily_survival_function(df = y.low,
                                            interval_name = interval.long,
                                            max_interval_val = max_interval_long,
                                            other_int_name = 'interval.short',
                                            other_int_fate = 'fate.short',
                                            interval_fate = 'fate.long')

y.med.daily.short <- daily_survival_function(df = y.med,
                                             interval_name = interval.short,
                                             max_interval_val = max_interval_short,
                                             other_int_name = 'interval.long',
                                             other_int_fate = 'fate.long',
                                             interval_fate = 'fate.short')

y.med.daily.long <- daily_survival_function(df = y.med,
                                            interval_name = interval.long,
                                            max_interval_val = max_interval_long,
                                            other_int_name = 'interval.short',
                                            other_int_fate = 'fate.short',
                                            interval_fate = 'fate.long')

y.high.daily.short <- daily_survival_function(df = y.high,
                                             interval_name = interval.short,
                                             max_interval_val = max_interval_short,
                                             other_int_name = 'interval.long',
                                             other_int_fate = 'fate.long',
                                             interval_fate = 'fate.short')

y.high.daily.long <- daily_survival_function(df = y.high,
                                            interval_name = interval.long,
                                            max_interval_val = max_interval_long,
                                            other_int_name = 'interval.short',
                                            other_int_fate = 'fate.short',
                                            interval_fate = 'fate.long')

# Export daily datasets ---------------------------------------------------------

write.csv(y.low.daily.short, here("data_outputs",
                          "03_simulated",
                          "02_analysis_ready",
                          "low_var_daily_short_data.csv"),
          row.names = F)

write.csv(y.low.daily.long, here("data_outputs",
                                 "03_simulated",
                                  "02_analysis_ready",
                                  "low_var_daily_long_data.csv"),
          row.names = F)

write.csv(y.med.daily.short, here("data_outputs",
                                  "03_simulated",
                          "02_analysis_ready",
                          "med_var_daily_short_data.csv"),
          row.names = F)

write.csv(y.med.daily.long, here("data_outputs",
                                 "03_simulated",
                                  "02_analysis_ready",
                                  "med_var_daily_long_data.csv"),
          row.names = F)

write.csv(y.high.daily.short, here("data_outputs",
                                   "03_simulated",
                           "02_analysis_ready",
                           "high_var_daily_short_data.csv"),
          row.names =F)

write.csv(y.high.daily.long, here("data_outputs",
                                  "03_simulated",
                                   "02_analysis_ready",
                                   "high_var_daily_long_data.csv"),
          row.names =F)

# Interval summarized data ------------------------------------------------

# Get interval summarised covariates --------------------------------------

interval_x_fun <- function(df, 
                           interval_name){
  
  df2 <- df %>%
    filter(day <= max_day) %>%
    distinct(Dataset, {{interval_name}}, day, x) %>%
    group_by(Dataset, {{interval_name}}) %>%
    summarise(x = mean(x)) %>%
    ungroup() %>%
    rename('interval' = {{interval_name}})%>%
    
  return(df2)
  
}

#low
y.low.short.x <- interval_x_fun(df = y.low,
                                interval_name = interval.short)

y.low.long.x <- interval_x_fun(df = y.low,
                               interval_name = interval.long)

#med
y.med.short.x <- interval_x_fun(df = y.med,
                                interval_name = interval.short)

y.med.long.x <- interval_x_fun(df = y.med,
                               interval_name = interval.long)
#high
y.high.short.x <- interval_x_fun(df = y.high,
                                interval_name = interval.short)

y.high.long.x <- interval_x_fun(df = y.high,
                               interval_name = interval.long)

# Get interval data summarised --------------------------------------------
  
interval_y_function <- function(df,
                                x.df){
  df2 <- df %>%
    filter(!is.na(fate)) %>%
    dplyr::select(Dataset, ID, interval, fate) %>%
    group_by(Dataset, ID, interval) %>%
    mutate(t = n()) %>%
    ungroup() %>%
    distinct(Dataset, ID, interval, fate,  t) %>%
    left_join(x.df, by = c('Dataset', 'interval'))
  
  return(df2)
  
}

#low
y.low.int.short <- interval_y_function(df = y.low.daily.short,
                    x.df = y.low.short.x)

y.low.int.long <- interval_y_function(df = y.low.daily.long,
                                       x.df = y.low.long.x)

#med
y.med.int.short <- interval_y_function(df = y.med.daily.short,
                                       x.df = y.med.short.x)

y.med.int.long <- interval_y_function(df = y.med.daily.long,
                                      x.df = y.med.long.x)

#high
y.high.int.short <- interval_y_function(df = y.high.daily.short,
                                       x.df = y.high.short.x)

y.high.int.long <- interval_y_function(df = y.high.daily.long,
                                      x.df = y.high.long.x)
# Export interval data ----------------------------------------------------

#low
write.csv(y.low.int.short, here("data_outputs",
                          "03_simulated",
                          "02_analysis_ready",
                          "low_var_short_interval_data.csv"),
          row.names =F)

write.csv(y.low.int.long, here("data_outputs",
                                "03_simulated",
                                "02_analysis_ready",
                                "low_var_long_interval_data.csv"),
          row.names =F)

#med
write.csv(y.med.int.short, here("data_outputs",
                                "03_simulated",
                                "02_analysis_ready",
                                "med_var_short_interval_data.csv"),
          row.names =F)

write.csv(y.med.int.long, here("data_outputs",
                               "03_simulated",
                               "02_analysis_ready",
                               "med_var_long_interval_data.csv"),
          row.names =F)

#high
write.csv(y.high.int.short, here("data_outputs",
                                "03_simulated",
                                "02_analysis_ready",
                                "high_var_short_interval_data.csv"),
          row.names =F)

write.csv(y.high.int.long, here("data_outputs",
                               "03_simulated",
                               "02_analysis_ready",
                               "high_var_long_interval_data.csv"),
          row.names =F)

# Full survey summarised data ---------------------------------------------

#for the total exposure model, need the above data to be summarised
#across the entire "lifetime" of the indiviudal, which is what
#these next steps do

#should this be... until the end of the survey interval as if we 
#went back only once, or should it be by individual?
#i think the first... yessss

total_y_function <- function(df){
  
  df2 <- df %>%
    filter(!is.na(fate)) %>%
    #for each dataset and id in each dataset
    group_by(Dataset, ID) %>%
    #get fate to be 0 if any are 0, otherwise 1
    mutate(fate = case_when(any(fate == 0) ~ 0,
                            TRUE ~ 1),
           #total number of "time" (equal to number of days)
           #setting to equal to the longest time period for all
           #such that we just went back for each individual after that
           #amount of time for all individuals 
           t= max_day,
           #x is mean of the value of x for that individual
           x = mean(x)) %>%
    ungroup() %>%
    #select dataset, id, fate, t, x
    distinct(Dataset, ID, fate, t, x)
  
  return(df2)
  
}

#i actually think these are the same regardless of 
#interval length??
low.var.tot <- total_y_function(y.low.daily.short)
med.var.tot <- total_y_function(y.med.daily.short)
high.var.tot <- total_y_function(y.high.daily.short)

# Export total exposure datasets ------------------------------------------

write.csv(low.var.tot, here("data_outputs",
                            "03_simulated",
                            "02_analysis_ready",
                            "low_var_total_data.csv"),
          row.names = F)

write.csv(med.var.tot, here("data_outputs",
                            "03_simulated",
                            "02_analysis_ready",
                            "med_var_total_data.csv"),
          row.names = F)

write.csv(high.var.tot, here("data_outputs",
                             "03_simulated",
                             "02_analysis_ready",
                             "high_var_total_data.csv"),
          row.names = F)

