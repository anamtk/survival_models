# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Data specifications -----------------------------------------------------


y_function <- function(n.ind, 
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
    rownames_to_column(var = "ID") %>%
    pivot_longer(V1:V10,
                 names_to = "interval",
                 values_to = "fate") 
  
  #which individuals are alive at the end of the intervals?
  #subset just those individuals in a dataframe
  alive <- df1 %>%
    group_by(ID) %>%
    filter(all(fate == 1)) %>%
    ungroup() 
  
  #which individuals are dead at some point in any
  #interval?
  #filter those out and then filter out any intervals
  #after the first one in which they are dead
  dead <- df1 %>%
    group_by(ID) %>%
    filter(any(fate == 0)) %>%
    mutate(first_0 = fate == 0 & !duplicated(fate == 0)) %>%
    # Find the first row with `double_B` in each group, filter out rows after it
    filter(row_number() <= min(which(fate == 0 & first_0 == TRUE))) %>%
    ungroup() %>%
    dplyr::select(-first_0)
  
  #combine those and turn back into a wide format matrix
  matrix2 <- alive %>%
    bind_rows(dead) %>%
    pivot_wider(names_from = "interval",
                values_from = "fate") %>%
      mutate(ID = as.numeric(ID)) %>%
      arrange(ID) %>%
    column_to_rownames(var = "ID") %>%
    as.matrix()
  
  return(matrix2)
  
}

# Getting y data prediction -----------------------------------------------

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
    rownames_to_column(var = "ID") %>%
    pivot_longer(V1:V10,
                 names_to = "interval",
                 values_to = "fate") %>%
    mutate(interval = str_sub(interval, start = 2, end = nchar(interval))) %>%
    mutate(interval = as.numeric(interval))
  
  #which individuals are alive at the end of the intervals?
  #subset just those individuals in a dataframe
  alive <- df1 %>%
    group_by(ID) %>%
    filter(all(fate == 1)) %>%
    ungroup() 
  
  #which individuals are dead at some point in any
  #interval?
  #filter those out and then filter out any intervals
  #after the first one in which they are dead
  dead <- df1 %>%
    group_by(ID) %>%
    filter(any(fate == 0)) %>%
    mutate(first_0 = fate == 0 & !duplicated(fate == 0)) %>%
    # Find the first row with `double_B` in each group, filter out rows after it
    filter(row_number() <= min(which(fate == 0 & first_0 == TRUE))) %>%
    ungroup() %>%
    dplyr::select(-first_0)
  
  #combine those and turn back into a wide format matrix
  matrix2 <- alive %>%
    bind_rows(dead) 
  
  return(matrix2)
  
}


# Get n.t for model3 ------------------------------------------------------

#filter data by dataset
#sort by 1's then >1's
#pull out the y's, t's, n.t's, and x's in that order
#make these into matrices


n.t_fun <- function(df, dataset){
  
  df1 <- df %>%
    filter(Dataset == dataset) %>%
    group_by(ID) %>%
    mutate(n.t = n()) %>%
    arrange(n.t) %>%
    ungroup() %>%
    distinct(ID, n.t)
  
  n.t <- as.vector(df1$n.t)
  
  return(n.t)
}

# Get X by dataset for model3 ---------------------------------------------

x_fun <- function(df, dataset){
  
  df1 <- df %>%
    filter(Dataset == dataset) %>%
    group_by(ID) %>%
    mutate(n.t = n()) %>%
    arrange(n.t) %>%
    ungroup() %>%
    dplyr::select(ID, interval, x) %>%
    pivot_wider(names_from = interval,
                values_from = x) %>%
    dplyr::select(-ID) %>%
    as.matrix()
  
  return(df1)
}

# Get t by dataset for model3 ---------------------------------------------


t_fun <- function(df, dataset){
  
  df1 <- df %>%
    filter(Dataset == dataset) %>%
    group_by(ID) %>%
    mutate(n.t = n()) %>%
    arrange(n.t) %>%
    ungroup() %>%
    dplyr::select(ID, interval, t) %>%
    pivot_wider(names_from = interval,
                values_from = t) %>%
    dplyr::select(-ID) %>%
    as.matrix()
  
  return(df1)
  
}

# Get y by datset for model3 ----------------------------------------------

y_fun <- function(df, dataset){
  
  y <- df %>%
    filter(Dataset == dataset) %>%
    group_by(ID) %>%
    mutate(n.t = n()) %>%
    arrange(n.t) %>%
    filter(interval == max(interval)) %>%
    ungroup() %>%
    dplyr::select(fate) %>%
    as_vector()
  
  return(y)
}


# Posterior dataframe from all models -------------------------------------

post_fun <- function(mod1, mod2, mod3){
  
  mod1_bs <- as.data.frame(mod1$quantiles) %>%
    rownames_to_column(var = "param") %>%
    dplyr::select(param, `50%`) %>%
    filter(param != "deviance") %>%
    mutate(parameter = str_sub(param, start = 1, end = 2)) %>%
    mutate(dataset = str_sub(param, start = 4, nchar(param)-1)) %>%
    dplyr::select(-param) %>%
    pivot_wider(names_from = parameter, 
                values_from = `50%`) %>%
    mutate(model = "Total Exposure")
  
  mod2_bs <- as.data.frame(mod2$quantiles) %>%
    rownames_to_column(var = "param") %>%
    dplyr::select(param, `50%`) %>%
    filter(param != "deviance") %>%
    mutate(parameter = str_sub(param, start = 1, end = 2)) %>%
    mutate(dataset = str_sub(param, start = 4, nchar(param)-1)) %>%
    dplyr::select(-param) %>%
    pivot_wider(names_from = parameter, 
                values_from = `50%`) %>%
    mutate(model = "Interval")
  
  mod3_bs <- as.data.frame(mod3$quantiles) %>%
    rownames_to_column(var = "param") %>%
    dplyr::select(param, `50%`) %>%
    filter(param != "deviance") %>%
    mutate(parameter = str_sub(param, start = 1, end = 2)) %>%
    mutate(dataset = str_sub(param, start = 4, nchar(param)-1)) %>%
    dplyr::select(-param) %>%
    pivot_wider(names_from = parameter, 
                values_from = `50%`) %>%
    mutate(model = "Custom")
  
  mod_results <- mod1_bs %>%
    bind_rows(mod2_bs, mod3_bs)
  
  return(mod_results)
  
}


# Percent post ------------------------------------------------------------

#how many of posteriors contain the true value?

post_perc_fun <- function(mod1, mod2, mod3, parm, beta){
  
  df1 <- as.data.frame(mod1$quantiles) %>%
    rownames_to_column(var = "param") %>%
    filter(param != "deviance") %>%
    mutate(parameter = str_sub(param, start = 1, end = 2)) %>%
    mutate(dataset = str_sub(param, start = 4, nchar(param)-1)) %>%
    dplyr::select(-param) %>%
    filter(parameter == beta) %>%
    mutate(category = case_when(parm < `2.5%` ~ "below",
                                parm > `97.5%` ~ "above",
                                (`2.5%` < parm & parm < `97.5%`) ~ "within")) %>%
    mutate(model = "Total Exposure")
  
  df2 <- as.data.frame(mod2$quantiles) %>%
    rownames_to_column(var = "param") %>%
    filter(param != "deviance") %>%
    mutate(parameter = str_sub(param, start = 1, end = 2)) %>%
    mutate(dataset = str_sub(param, start = 4, nchar(param)-1)) %>%
    dplyr::select(-param) %>%
    filter(parameter == beta) %>%
    mutate(category = case_when(parm < `2.5%` ~ "below",
                                parm > `97.5%` ~ "above",
                                (`2.5%` < parm & parm < `97.5%`) ~ "within")) %>%
    mutate(model = "Interval")
  
  df3 <- as.data.frame(mod3$quantiles) %>%
    rownames_to_column(var = "param") %>%
    filter(param != "deviance") %>%
    mutate(parameter = str_sub(param, start = 1, end = 2)) %>%
    mutate(dataset = str_sub(param, start = 4, nchar(param)-1)) %>%
    dplyr::select(-param) %>%
    filter(parameter == beta) %>%
    mutate(category = case_when(parm < `2.5%` ~ "below",
                                parm > `97.5%` ~ "above",
                                (`2.5%` < parm & parm < `97.5%`) ~ "within")) %>%
    mutate(model = "Custom")
  
  all <- df1 %>%
    bind_rows(df2, df3)
  
  sum <- all %>%
    group_by(model, category) %>%
    tally() %>%
    ungroup() %>%
    pivot_wider(names_from = category,
                values_from = n) %>%
    rowwise() %>%
    replace_na(list(below = 0, within = 0, above = 0)) %>%
    mutate(percent = within/(within + below + above))
  
  return(sum)
  
}

