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

