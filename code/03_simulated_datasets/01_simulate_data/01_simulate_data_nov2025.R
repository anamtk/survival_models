
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
b0 <- 5.5
b1 <- 1.25

#number of simulated datasets
n.datasets <- 100

#UPDATE:::
#simulate x data for every day of every interval
#in every dataset and KEEP them all

#Then generate y data for those x's as we've done - then, 
#"fill" 0's for intervals in which the individual "died"
#since this is the way the surveys would have been performed
#so for the last interval, we get all days in that interval, 
#and all are set to 0 - since we don't know which day in
#the total survey interval an individual died. 

#then combine:
#1. daily x data for each individual in each dataset
#in each day with interval IDs recoginzed
#2. interval level survival data with dataset, ID, and 
#interval ID info

sim_x_function <- function(n.t,
                           var_level,
                           b0, 
                           b1){
  
  #creatng the x so that it gradually changes 
  #throughout the series
  days <- as.data.frame(1:n.t)
  
  names(days) <- "day"
  
  days2 <- days %>%
    rowwise() %>%
    #this generally goes from -sd to +sd through 45 days
    mutate(x1 = rnorm(n(), mean(0.01*day - 0.1), sd =1))
  
  #get x for each interval (day)
  #x1 <- rnorm(n.t, mean = 0, sd = 1)
  #scale that 
  x2 <- scale(days2$x1)
  #multiply by the level of variation we want
  #in our dataset (low 0.01, med 0.1, high 1)
  x.vec <- x2*var_level
  
  #divide it into intervals of ~3.5 days each
  interval <- tibble(interval = sort(rep(1:15, 3))) %>%
    mutate(day = 1:n())
  
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
    #filter IDs where all intervals == 1
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
    dplyr::select(ID, day, Dataset, interval, fate)
  
  ind.df2 <- ind.df %>%
    dplyr::select(-fate) %>%
    left_join(df2, by = c("ID", "day", "Dataset",
                          "interval"))%>%
    group_by(ID, interval, Dataset) %>%
    mutate(fate = case_when(any(fate == 0) ~ 0,
                            TRUE ~ fate)) %>%
    ungroup() %>%
    filter(!is.na(fate))
  
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
                                 var_level = 0.01) #based on this probability
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
                                 var_level = 0.1) #based on this probability
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
                                 var_level = 1) #based on this probability
                } )


x.high.df <- bind_rows(x.high, .id = "Dataset") 

set.seed(1)
y.high <- sim_y_function(n.ind = n.indiv,
                        n.t = n.int,
                        n.data = n.datasets,
                        x.df = x.high.df)

# Explore datasets --------------------------------------------------------
theme_set(theme_bw())
c <- y.high %>%
  ggplot(aes(x = day, y= x, group = Dataset)) +
  geom_line() +
  geom_vline(aes(xintercept = 3), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 6), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 9), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 12), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 15), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 18), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 21), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 24), linetype = 2, color = "blue") +
  xlim(0, 36)

b <- y.med %>%
  ggplot(aes(x = day, y= x, group = Dataset)) +
  geom_line()+
  geom_vline(aes(xintercept = 3), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 6), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 9), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 12), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 15), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 18), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 21), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 24), linetype = 2, color = "blue") +
  xlim(0, 36)

a <- y.low %>%
  ggplot(aes(x = day, y= x, group = Dataset)) +
  geom_line()+
  geom_vline(aes(xintercept = 3), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 6), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 9), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 12), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 15), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 18), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 21), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 24), linetype = 2, color = "blue") +
  xlim(0, 36)

#a+ b+ c

d <- y.low %>%
  ggplot(aes(x = day, y = ps, group = Dataset)) +
  geom_line()+
  geom_vline(aes(xintercept = 3), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 6), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 9), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 12), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 15), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 18), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 21), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 24), linetype = 2, color = "blue") +
  xlim(0, 36)

e <- y.med %>%
  ggplot(aes(x = day, y = ps, group = Dataset)) +
  geom_line()+
  geom_vline(aes(xintercept = 3), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 6), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 9), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 12), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 15), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 18), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 21), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 24), linetype = 2, color = "blue") +
  xlim(0, 36)

f <- y.high %>%
  ggplot(aes(x = day, y = ps, group = Dataset)) +
  geom_line()+
  geom_vline(aes(xintercept = 3), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 6), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 9), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 12), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 15), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 18), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 21), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 24), linetype = 2, color = "blue") +
  xlim(0, 36)

(a+ b+ c)/
(d + e +f)

a+b+c
a

# Truncate so that a reasonable number in last interval -------------------

#for the high var, the max survival is to:
(max_day <- y.high %>%
  filter(!is.na(fate)) %>%
  filter(day == max(day)) %>%
   distinct(day) %>%
   dplyr::select(day) %>%
   as_vector())

#just setting below that so that we have a 
#fair number that make it to the last interval
max_day <- 23

#get that max interval
(max_interval <- y.high %>%
  filter(!is.na(fate)) %>%
  filter(day == max_day) %>%
  distinct(interval) %>%
  dplyr::select(interval) %>%
  as_vector() )

#get data truncated to that point for all variation levels
y.low.daily <- y.low %>%
  filter(interval <= max_interval) %>%
  filter(!is.na(fate))

y.med.daily <- y.med %>%
  filter(interval <= max_interval) %>%
  filter(!is.na(fate))

y.high.daily <- y.high %>%
  filter(interval <= max_interval) %>%
  filter(!is.na(fate))

#make sure some in the high variability
#survive to the last interval
#try to get the class imbalance??
#aiming for ~20% dying if basing on the 
#woodpecker dataset
#that's 600 in 0 class
#out of 3000 individuals,
#or ~6 individuals of every set of 30
y.high.daily %>%
  group_by(Dataset, ID) %>%
  filter(interval == max(interval)) %>%
  distinct(Dataset, ID, fate) %>%
  group_by(Dataset, fate) %>%
  tally() %>%
  mutate(fate = as.factor(fate)) %>%
  ggplot(aes(x = fate, y = n, fill = fate))+
  geom_boxplot()
#looks pretty close

y.high.daily %>%
  group_by(Dataset, ID) %>%
  filter(interval == max(interval)) %>%
  distinct(Dataset, ID, fate) %>%
  group_by(fate) %>%
  tally()


# Redistribute "live" individuals -----------------------------------------

#right now, the only "live" individuals made it 
#to the last interval, so we need to redistribute those 
#live intervals so that there are 1-0 data at every
#set of interval values. I used the distribution of 
#values across the white-headed woodpecker dataset for 
#the number of intervals to define this, but it generally
#peaks at 2-3 resurveys and tapers off, so I'm going to take
#all those indivdiausl that "lived" to the last interval and 
#just pretend we only surveyed them for one of 1-7 resurveyes
#based on that distribution

reassign_live_fun <- function(df, dataset){
  
  #create an "alive" dataset 
  alive <- df %>%
    filter(Dataset == dataset) %>%
    group_by(ID, Dataset) %>%
    filter(all(fate == 1)) %>%
    ungroup() 
  
  #just live IDs 
  ids <- alive %>%
    distinct(ID, Dataset)
  
  #get the total number so we can normalize 
  #the resampled number 
  total <- nrow(ids)
  
  #intervals to re-assign individuals to
  intervals <- 1:(max_interval-1)
  
  #this is assuming that there are only 7 intervals
  #which i think is fine for now, but would need to update
  #if we added more intervals down the road
  #gives them a humped distribution where more are surveyed
  #for 2-3 times, which is what was the distribution
  #in the woodpecker dataset
  reassign_alive <- c(rep(intervals[1], floor(0.13*total)),
                      rep(intervals[2], floor(0.2*total)),
                      rep(intervals[3], floor(0.2*total)),
                      rep(intervals[4], floor(0.17*total)),
                      rep(intervals[5], floor(0.13*total)),
                      rep(intervals[6], floor(0.1*total)),
                      rep(intervals[7], floor(0.1*total)))
  
  pad_value <- 8
  
  #in case that vector isn't the same length as the number
  #of live individuals, pad with interval 8 at the end
  #to the length - those won't be re-assigned
  reassign_alive <- rep_len(c(reassign_alive, pad_value),
                            length.out = total)
  
  #give each ID a new max interval
  ids2 <- ids %>%
    bind_cols(reassign_alive) %>%
    rename("new_max_interval" = "...3")
  
  #bind with the "live" dataset and filter that dataset
  #to be just those intervals
  alive1 <- alive %>%
    left_join(ids2, by = c("ID", "Dataset")) %>%
    group_by(ID, Dataset) %>%
    filter(interval <= new_max_interval)
  
  #create an "dead" dataset 
  dead <- y.low.daily %>%
    filter(Dataset == dataset) %>%
    group_by(ID, Dataset) %>%
    filter(any(fate == 0)) %>%
    ungroup() 
  
  all_df <- alive1 %>%
    bind_rows(dead)
  
  return(all_df)
}

#low variation reassignment of live individuals:
y.low.daily.list2 <- lapply(1:100, #how many datasets
                 function(i)
                 {
                   reassign_live_fun(df = y.low.daily,
                                     dataset = i) #based on this probability
                 } )

y.low.daily2 <- bind_rows(y.low.daily.list2)

#med variation reassignment
y.med.daily.list2 <- lapply(1:100, #how many datasets
                            function(i)
                            {
                              reassign_live_fun(df = y.med.daily,
                                                dataset = i) #based on this probability
                            } )

y.med.daily2 <- bind_rows(y.med.daily.list2)

#high
y.high.daily.list2 <- lapply(1:100, #how many datasets
                            function(i)
                            {
                              reassign_live_fun(df = y.high.daily,
                                                dataset = i) #based on this probability
                            } )

y.high.daily2 <- bind_rows(y.high.daily.list2)

# Export daily datasets ---------------------------------------------------------

write.csv(y.low.daily2, here("data_outputs",
                          "03_simulated",
                          "02_analysis_ready",
                          "low_var_daily_data.csv"),
          row.names = F)

write.csv(y.med.daily2, here("data_outputs",
                          "03_simulated",
                          "02_analysis_ready",
                          "med_var_daily_data.csv"),
          row.names = F)

write.csv(y.high.daily2, here("data_outputs",
                           "03_simulated",
                           "02_analysis_ready",
                           "high_var_daily_data.csv"),
          row.names = F)

# Interval summarized data ------------------------------------------------

# Get interval summarised covariates --------------------------------------

y.low.x <- y.low.daily2 %>%
  distinct(Dataset, interval, day, x) %>%
  group_by(Dataset, interval) %>%
  summarise(x = mean(x)) %>%
  ungroup()

y.med.x <- y.med.daily2 %>%
  distinct(Dataset, interval, day, x) %>%
  group_by(Dataset, interval) %>%
  summarise(x = mean(x)) %>%
  ungroup()

y.high.x <- y.high.daily2 %>%
  distinct(Dataset, interval, day, x) %>%
  group_by(Dataset, interval) %>%
  summarise(x = mean(x)) %>%
  ungroup()

y.low.int <- y.low.daily2 %>%
  filter(!is.na(fate)) %>%
  dplyr::select(Dataset, ID, interval, fate) %>%
  group_by(Dataset, ID, interval) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ fate)) %>%
  ungroup() %>%
  group_by(Dataset, ID, interval) %>%
  mutate(t = n()) %>%
  ungroup() %>%
  distinct(Dataset, ID, interval, fate, t) %>%
  left_join(y.low.x, by = c('Dataset', 'interval'))


y.med.int <- y.med.daily2 %>%
  filter(!is.na(fate)) %>%
  dplyr::select(Dataset, ID, interval, fate) %>%
  group_by(Dataset, ID, interval) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ fate)) %>%
  ungroup() %>%
  group_by(Dataset, ID, interval) %>%
  mutate(t = n()) %>%
  ungroup() %>%
  distinct(Dataset, ID, interval, fate, t) %>%
  left_join(y.med.x, by = c('Dataset', 'interval'))

y.high.int <- y.high.daily2 %>%
  filter(!is.na(fate)) %>%
  dplyr::select(Dataset, ID, interval, fate) %>%
  group_by(Dataset, ID, interval) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ fate)) %>%
  ungroup() %>%
  group_by(Dataset, ID, interval) %>%
  mutate(t = n()) %>%
  ungroup() %>%
  distinct(Dataset, ID, interval, fate, t) %>%
  left_join(y.high.x, by = c('Dataset', 'interval'))
  

# Export interval data ----------------------------------------------------

write.csv(y.low.int, here("data_outputs",
                          "03_simulated",
                          "02_analysis_ready",
                          "low_var_interval_data.csv"),
          row.names = F)

write.csv(y.med.int, here("data_outputs",
                          "03_simulated",
                          "02_analysis_ready",
                          "med_var_interval_data.csv"),
          row.names = F)

write.csv(y.high.int, here("data_outputs",
                           "03_simulated",
                           "02_analysis_ready",
                           "high_var_interval_data.csv"),
          row.names = F)


# Full survey summarised data ---------------------------------------------

#for the total exposure model, need the above data to be summarised
#across the entire "lifetime" of the indiviudal, which is what
#these next steps do


low.var.tot <- y.low.daily2 %>%
  filter(!is.na(fate)) %>%
  #for each dataset and id in each dataset
  group_by(Dataset, ID) %>%
  #get fate to be 0 if any are 0, otherwise 1
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         #total number of "time" (equal to number of days)
         t= n(),
         #x is mean of the value of x for that individual
         x = mean(x)) %>%
  ungroup() %>%
  #select dataset, id, fate, t, x
  distinct(Dataset, ID, fate, t, x)

#do the same process for the medium variability dataset
med.var.tot <- y.med.daily2 %>%
  group_by(Dataset, ID) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         t= n(),
         x = mean(x)) %>%
  ungroup() %>%
  distinct(Dataset, ID, fate, t, x)

#do the same for the high variability dataset
high.var.tot <- y.high.daily2 %>%
  group_by(Dataset, ID) %>%
  mutate(fate = case_when(any(fate == 0) ~ 0,
                          TRUE ~ 1),
         t= n(),
         x = mean(x)) %>%
  ungroup() %>%
  distinct(Dataset, ID, fate, t, x)

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

