set.seed(1)
#number of individuals in dataset
n.indiv <- 300
#number of times each individual was smpaled, with max 4
n.t <- sample(4, 300, replace = T)

#get a metadata file of the total number of surveys/individual
meta <- as.data.frame(cbind(ID = 1:300,
                    n.t = n.t))


#get full visit dataset by repeating id by the n.t for the ID
visits <- as.data.frame(rep(meta$ID, meta$n.t)) %>%
  rename("ID" = 'rep(meta$ID, meta$n.t)')%>%
  group_by(ID) %>%
  #get a 1:n interval number per ID
  mutate(int = 1:n()) %>%
  ungroup() 

#mean of 9 day interval from random samples
t <- round(rnorm(nrow(visits), mean = 9, sd = 3))
t <- as.data.frame(t)

visits2 <- visits %>%
  bind_cols(t) %>%
  group_by(ID) %>%
  mutate(end = cumsum(t)) %>%
  mutate(start = lag(end)) %>%
  mutate(start = case_when(is.na(start) ~ 1,
                           TRUE ~ start)) %>%
  ungroup()
  
#what is the maximum time an individual was surveyed?
visits2 %>%
  group_by(ID) %>%
  summarise(sum = sum(t)) %>%
  summarise(max = max(sum)) #49 days is max in this dataset right now
  
#does variable need to vary with time? or does it just need
#to vary in the season in low, med, high levels? HMM...
time <- 1:49
x1 <- time*5 + 20

plot(x1 ~ time)

low <- x1 + runif(length(x1), min = 0, max = 3)
low <- as.data.frame(low) %>%
  mutate(day = 1:n())

med <- x1 + runif(length(x1), min = 0, max = 8)
med <- as.data.frame(med) %>%
  mutate(day = 1:n())


high <- x1 + runif(length(x1), min = 0, max = 15)
high <- as.data.frame(high) %>%
  mutate(day = 1:n())


covs <- low %>%
  left_join(med, by = "day") %>%
  left_join(high, by = "day")

# Function to filter daily covariate dataset by ID and time interval
cov_int_mean <- function(start_date, end_date, level){
  
  dat <- covs %>%
    filter(day >= start_date & day < end_date)
  
  var <- dat %>%
    dplyr::select({{level}}) %>%
    summarise(mean = mean({{level}})) %>%
    as_vector()
  
  return(var)
  
}


cov_int_mean(start_date = 3,
             end_date = 7, 
             level = high)

visits2$low <- rep(NA, length(visits2$ID))

for(i in 1:length(visits2$ID)) {
  visits2$low[i] <- cov_int_mean(start_date = visits2$start[i], 
                                 end_date = visits2$end[i],
                                 level = low) 
}

visits2$med <- rep(NA, nrow(visits2))

for(i in 1:nrow(visits2)){
  visits2$med[i] <- cov_int_mean(start_date = visits2$start[i], 
                                 end_date = visits2$end[i],
                                 level = med) 
}

visits2$high <- rep(NA, nrow(visits2))

for(i in 1:nrow(visits2)){
  visits2$high[i] <- cov_int_mean(start_date = visits2$start[i], 
                                 end_date = visits2$end[i],
                                 level = high) 
}

hist(visits2$low) 
hist(visits2$med) 
hist(visits2$high) 

visits2 <- visits2 %>%
  mutate(low_st = scale(low),
         med_st = scale(med),
         high_st = scale(high))

beta <- 0.5

z = beta*(visits2$low_st)        # linear combination with a bias
pr = exp(z)/(1+exp(z))         # pass through an inv-logit function
visits2$y = rbinom(nrow(visits2),1,pr) 

