# Predictive accuracy figures
# May 4, 2023
# Ana Miller-ter Kuile

# this is a visualization of the predictive accuracy for each
# of the three models for the nest dataset

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse",
                  "reshape2", 
                  "pROC",
                  "patchwork")


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

# Predictive accuracy -----------------------------------------------------

# Model 1 -----------------------------------------------------------------

mu_p1 <- as.data.frame(mod1_GOF$mean$p) %>%
  rename("P" = 'mod1_GOF$mean$p')

y_acc1 <- y1 %>%
  bind_cols(mu_p1) %>%
  mutate(Fate_class = as.factor(Fate_class)) %>%
  mutate(type = case_when((Fate_class == 0 & P >= .5) ~ "mis-0",
                          (Fate_class == 0 & P < .5) ~ "match-0",
                          (Fate_class == 1 & P >= 0.5) ~ "match-1",
                          (Fate_class == 1 & P < 0.5) ~ "mis-1"))

y_acc1 %>%
  group_by(type) %>%
  tally()

#accuary
#0s:
zeros1 <- round(79/(79+13), digits = 2)*100
#1s:
ones1 <- round(225/(225+3), digits = 2)*100

(mod1_acc_plot <- ggplot(y_acc1, aes(x = Fate_class, y = P)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_boxplot(fill = '#8dd3c7') +
    labs(x = "Observed fate",
         y = "Predicted survival probability",
         title = "A. Total exposure model") +
    annotate(geom = "text", 
             x = 0.75, y = 0.45,
             label = paste(zeros1, "%", sep = "")) +
    annotate(geom = "text", 
             x = 2.25, y = 0.55,
             label = paste(ones1, "%", sep = "")) +
    ylim(0, 1) +
    theme(axis.title.x = element_blank()))

y_acc1 %>%
  group_by(Fate_class) %>%
  summarise(meanp = mean(P),
            sdp = sd(P),
            total = n(),
            sep = sdp/sqrt(total))

# Model 2 -----------------------------------------------------------------


mu.p2 <- as.data.frame(mod2_GOF$mean$p.intkeep) %>%
  rename("P" = 'mod2_GOF$mean$p.intkeep')

y_acc2 <- y2 %>%
  bind_cols(mu.p2) %>%
  mutate(Fate_class = as.factor(Fate_class)) %>%
  mutate(type = case_when((Fate_class == 0 & P >= .5) ~ "mis-0",
                          (Fate_class == 0 & P < .5) ~ "match-0",
                          (Fate_class == 1 & P >= 0.5) ~ "match-1",
                          (Fate_class == 1 & P < 0.5) ~ "mis-1"))


y_acc2 %>%
  group_by(type) %>%
  tally()

#1s:
ones2 <- round(228/228, digits = 2)*100
#0s:
zeros2 <- round(3/(3+89), digits = 2)*100

(mod2_acc_plot <- ggplot(y_acc2, aes(x = Fate_class, y = P)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_boxplot(fill = "#bebada") +
    labs(x = "Observed fate",
         y = "Predicted survival probability",
         title = "B. Interval data model") +
    annotate(geom = "text", 
             x = 0.75, y = 0.45,
             label = paste(zeros2, "%", sep = "")) +
    annotate(geom = "text", 
             x = 2.25, y = 0.55,
             label = paste(ones2, "%", sep = "")) +
    ylim(0, 1) +
    theme(axis.title.y = element_blank()))

y_acc2 %>%
  group_by(Fate_class) %>%
  summarise(meanp = mean(P),
            sdp = sd(P),
            total = n(),
            sep = sdp/sqrt(total))


# Model 3 -----------------------------------------------------------------


mu_p3.1 <- as.data.frame(mod3_GOF$mean$yrep_1) %>%
  rename("P" = 'mod3_GOF$mean$yrep_1')

mu_p3.2 <- as.data.frame(mod3_GOF$mean$yrep_2) %>%
  rename("P" = 'mod3_GOF$mean$yrep_2') %>%
  filter(!is.na(P))

times <- as.data.frame(data3$n.t) %>%
  rename("intervals" = "data3$n.t")

mu_p3 <- mu_p3.1 %>%
  bind_rows(mu_p3.2)

y_acc3 <- y3 %>%
  bind_cols(mu_p3, times) %>%
  mutate(Fate_class = as.factor(Fate_class)) %>%
  mutate(type = case_when((Fate_class == 0 & P >= .5) ~ "mis-0",
                          (Fate_class == 0 & P < .5) ~ "match-0",
                          (Fate_class == 1 & P >= 0.5) ~ "match-1",
                          (Fate_class == 1 & P < 0.5) ~ "mis-1"))

y_acc3 %>%
  group_by(type) %>%
  tally()

#accuary
#0s:
zeros3 <- round(80/(80+12), digits = 2)*100#87%
#1s:
ones3 <- round(224/(224+4), digits = 2)*100 #98%

(mod3_acc_plot <- ggplot(y_acc3, aes(x = Fate_class, y = P)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_boxplot(fill = "#fdb462") +
    labs(x = "Observed fate",
         y = "Predicted survival probability",
         title = "C. Custom model")  +
    annotate(geom = "text", 
             x = 0.75, y = 0.45,
             label = paste(zeros3, "%", sep = "")) +
    annotate(geom = "text", 
             x = 2.25, y = 0.55,
             label = paste(ones3, "%", sep = "")) +
    ylim(0, 1) +
    theme(axis.title = element_blank()))

y_acc3 %>%
  group_by(Fate_class) %>%
  summarise(meanp = mean(P),
            sdp = sd(P),
            total = n(),
            sep = sdp/sqrt(total))



# Combine plots -----------------------------------------------------------

mod1_acc_plot + mod2_acc_plot + mod3_acc_plot

ggsave(filename = here("pictures",
                       "whwo_accuracy.pdf"),
       height = 4, 
       width = 8,
       units = "in")




