# Posterior predictive checks
# May 4, 2023
# Ana Miller-ter Kuile

#these are AUC distribution plots for the nest project

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "reshape2", "patchwork",
                  "pROC")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

source(here("code",
            "00_functions",
            "GOF_functions.R"))

source(here("code",
            "00_functions",
            "plot_functions.R"))


# Load GOF model runs -----------------------------------------------------

mod1_GOF <- readRDS(here('monsoon',
                         "01_whwonests",
                         "model1",
                         "outputs",
                         "model1_GOF.RDS"))


mod2_GOF <- readRDS(here('monsoon',
                         "01_whwonests",
                         "model2",
                         "outputs",
                         "model2_GOF.RDS"))

mod3_GOF <- readRDS(here('monsoon',
                         "01_whwonests",
                         "model3",
                         "outputs",
                         "model3_GOF.RDS"))

# Load data ---------------------------------------------------------------

#and we also need our original y data
data1 <- readRDS(here("data_outputs",
                      "01_whwonests",
                      "03_JAGS_input_data",
                      "mod1_JAGS_input_data.RDS"))

#and we also need our original y data
data2 <- readRDS(here("data_outputs",
                      "01_whwonests",
                      "03_JAGS_input_data",
                      "mod2_JAGS_input_data.RDS"))

data3 <- readRDS(here("data_outputs",
                      "01_whwonests",
                      "03_JAGS_input_data",
                      "mod3_JAGS_input_data.RDS"))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y1 <- as.data.frame(data1$y) %>%
  rename("Fate_class" = "data1$y") %>%
  mutate(Nest_ID = 1:n(),
         type = "Observed") 

#we need to extract our observed data from our dataframe
y2 <- as.data.frame(data2$y) %>%
  mutate(Nest_ID = 1:n()) %>%
  pivot_longer(1:15,
               names_to = "Interval",
               values_to = "Fate_class") %>%
  filter(!is.na(Fate_class)) %>%
  group_by(Nest_ID) %>%
  filter(Interval == max(Interval, na.rm = T)) %>%
  ungroup() %>%
  mutate(type = "Observed") 

#we need to extract our observed data from our dataframe
y3 <- as.data.frame(data3$y) %>%
  rename("Fate_class" = "data3$y") %>%
  mutate(Nest_ID = 1:n(),
         type = "Observed") 


# AUC ---------------------------------------------------------------------


# Model1 ------------------------------------------------------------------


resp1 <- as.vector(y1$Fate_class)

iteration.num1 <- length(mod1_GOF$sims.list$p[,1])

mod1_AUC <- rep(NA, iteration.num1)

for(i in 1:iteration.num1){
  mod1_AUC[i] <- AUC_JAGS(mod1_GOF, 
                          iteration.num = i, 
                          resp = resp1)
}

mean1 <- as.data.frame(mod1_AUC) %>%
  summarise(mean = mean(mod1_AUC)) %>%
  as_vector()

(mod1_AUC_plot <- as.data.frame(mod1_AUC) %>%
    ggplot() +
    geom_density(aes(x = mod1_AUC), fill = "#8dd3c7", color = "black") +
    geom_vline(xintercept = mean1, linetype = 2) +
    labs(x = "AUC",
         y = "Count",
         title = "A. Total exposure model") +
    xlim(0.4, 1) +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank()))


# Model2 ------------------------------------------------------------------


resp2 <- as.vector(y2$Fate_class)

iteration.num2 <- length(mod2_GOF$sims.list$p.intkeep[,1])

AUC_JAGS4(mod2_GOF, 
          iteration.num = 11, 
          resp = resp2)

mod2_AUC <- rep(NA, iteration.num2)

for(i in 1:iteration.num2){
  mod2_AUC[i] <- AUC_JAGS4(mod2_GOF, 
                           iteration.num = i, 
                           resp = resp2)
}

mean2 <- as.data.frame(mod2_AUC) %>%
  summarise(mean = mean(mod2_AUC)) %>%
  as_vector()

# (mod2_AUC_plot <- as.data.frame(mod2_AUC) %>%
#     ggplot() +
#     geom_histogram(aes(x = mod2_AUC)) +
#     geom_vline(xintercept = mean2, linetype = 2) +
#     labs(x = "AUC",
#          y = "Count",
#          title = "B. Interval data model"))

#THIS IS for all - would need to track pint though for all intervals
t <- mod2_GOF$sims.list$p.int

layers <- dim(t)[[3]]

dfs <- lapply(1:layers,
              function(x){
                return(as.data.frame(t[,,x]))
              } )

dfs1 <- dfs %>%
  map(~mutate(., iteration = 1:n()))

full_df <- bind_rows(dfs1, .id = "interval") %>%
  pivot_longer(cols = 2:(last_col()-1),
               values_to = "p",
               names_to = "Nest_ID") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, length(Nest_ID))) %>%
  unite(col = "ID_interval",
        c("Nest_ID", "interval"),
        sep = "_") %>%
  filter(!is.na(p))


resp2.1 <- as.data.frame(data2$y) %>%
  mutate(Nest_ID = 1:n()) %>%
  pivot_longer(cols = 1:(last_col()-1),
               values_to = "resp",
               names_to = "interval") %>%
  filter(!is.na(resp)) %>%
  unite(col = "ID_interval",
        c("Nest_ID", "interval"),
        sep = "_")


AUC_JAGS2(df = full_df,
          iteration.num = 3,
          resp = resp2.1$resp)

iteration.num2.1 <- length(unique(full_df$iteration))

mod2_AUC2 <- rep(NA, iteration.num2.1)

for(i in 1:iteration.num2.1){
  mod2_AUC2[i] <- AUC_JAGS2(df = full_df,
                            iteration.num = i,
                            resp = resp2.1$resp)
}

mean2.1 <- as.data.frame(mod2_AUC2) %>%
  summarise(mean = mean(mod2_AUC2)) %>%
  as_vector()

# (mod2_AUC_plotallints <- as.data.frame(mod2_AUC2) %>%
#     ggplot() +
#     geom_histogram(aes(x = mod2_AUC2)) +
#     geom_vline(xintercept = mean2.1, linetype = 2) +
#     labs(x = "AUC",
#          y = "Count",
#          title = "B. Interval data model"))

mod2auc <- as.data.frame(mod2_AUC) %>%
  rename("AUC" = "mod2_AUC") %>%
  mutate(type = "Last interval")
mod2auc2 <- as.data.frame(mod2_AUC2) %>%
  rename("AUC" = "mod2_AUC2") %>%
  mutate(type = "All intervals") %>%
  bind_rows(mod2auc)

(mod2_AUC_plotall <- mod2auc2 %>%
    ggplot() +
    geom_density(aes(x = AUC, fill = type), 
                   color = "black") +
    scale_fill_manual(values = c("#E1DFEE","#bebada")) +
    geom_vline(xintercept = mean2.1, linetype = 2) +
    geom_vline(xintercept = mean2, linetype = 2) +
    labs(x = "AUC",
         y = "Count",
         title = "B. Interval data model")+
    xlim(0.4, 1) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()))

# Model 3 -----------------------------------------------------------------

resp3 <- as.vector(data3$y)

iteration.num3 <- length(mod3_GOF$sims.list$p1[,1])

pred3.1 <- as.data.frame(t(mod3_GOF$sims.list$p1))
pred3.2 <- as.data.frame(t(mod3_GOF$sims.list$p2))
pred3.2 <- pred3.2[25:nrow(pred3.2),]

pred3 <- rbind(pred3.1, pred3.2)

mod3_AUC <- rep(NA, iteration.num3)

for(i in 1:iteration.num3){
  mod3_AUC[i] <- AUC_JAGS3(pred3, 
                           iteration.num = i, 
                           resp = resp3)
}

mean3 <- as.data.frame(mod3_AUC) %>%
  summarise(mean = mean(mod3_AUC)) %>%
  as_vector()

(mod3_AUC_plot <- as.data.frame(mod3_AUC) %>%
    ggplot() +
    geom_density(aes(x = mod3_AUC), fill = '#fdb462', color = "black") +
    geom_vline(xintercept = mean3, linetype = 2) +
    labs(x = "AUC",
         y = "Count",
         title = "C. Custom model") +
    xlim(0.4, 1) +
    theme(axis.title.y = element_blank()))


# Combine -----------------------------------------------------------------

mod1_AUC_plot + mod2_AUC_plotall + mod3_AUC_plot +
  plot_layout(ncol = 1)


