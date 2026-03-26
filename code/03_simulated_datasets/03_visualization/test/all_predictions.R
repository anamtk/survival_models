#Pull out median effects for each model
#Ana Miller-ter Kuile
#June 12, 2023  

#this script pulls out the posterior summaries for each of the models
#and looks at how well each model predicts the real values

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("code",
            "00_functions",
            "simulation_functions.R"))

theme_set(theme_bw())

b0_1 <- 2.1
b0_2 <- 3.5
b1 <- 3.5


# Load datasets -----------------------------------------------------------

file.names <- list.files(here('data_outputs',
                              '03_simulated',
                              '04_posterior_summaries',
                              'test'),
                         pattern = "summary_*.RDS")

t <- list.files(here('data_outputs',
                '03_simulated',
                '04_posterior_summaries',
                'test'),
           pattern = "summary.RDS",
           recursive = TRUE,
           full.names = TRUE) %>%
  map(readRDS)

names(t) <- file.names

stats_fun <- function(model){
  
  stats <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    mutate(type = case_when(str_detect(parm, "b1") ~ "b1",
                            str_detect(parm, "b0") ~ "b0",
                            TRUE ~ NA_character_))

  return(stats)
}

stats_list <- lapply(t, stats_fun)

stats_df <- bind_rows(stats_list, 
                      .id = "model") %>%
  mutate(model = str_sub(model, 1, (nchar(model)-12))) %>%
  filter(!is.na(type)) %>%
  dplyr::select(type,
              model,
              `2.5%`,
              `50%`,
              `97.5%`) %>%
  mutate(prediction = case_when((type == 'b1' & 
                                   `2.5%` < b1 &
                                   `97.5%` > b1) ~ "correct",
                                (type == "b1" &
                                   `2.5%` < b1 &
                                   `97.5%` < b1) ~ "incorrect",
                                (type == "b1" &
                                   `2.5%` > b1 &
                                   `97.5%` > b1) ~ "incorrect",
                                (str_detect(model, 'low1') &
                                   type == 'b0' & 
                                   `2.5%` < b0_1 &
                                   `97.5%` > b0_1) ~ "correct",
                                (str_detect(model, 'low1') &
                                   type == "b0" &
                                   `2.5%` < b0_1 &
                                   `97.5%` < b0_1) ~ "incorrect",
                                (str_detect(model, 'low1') &
                                   type == "b0" &
                                   `2.5%` > b0_1 &
                                   `97.5%` > b0_1) ~ "incorrect",
                                (str_detect(model, 'low2') &
                                   type == 'b0' & 
                                   `2.5%` < b0_2 &
                                   `97.5%` > b0_2) ~ "correct",
                                (str_detect(model, 'low2') &
                                   type == "b0" &
                                   `2.5%` < b0_2 &
                                   `97.5%` < b0_2) ~ "incorrect",
                                (str_detect(model, 'low1') &
                                   type == "b0" &
                                   `2.5%` > b0_2 &
                                   `97.5%` > b0_2) ~ "incorrect",
                                TRUE ~ NA_character_)) %>% 
  mutate(mod = str_sub(model, 1, 4),
         dataset = str_sub(model, nchar(model), nchar(model))) 
 
b01 <- stats_df %>%
  filter(type == "b0") %>%
  filter(dataset == "1") %>%
  group_by(prediction, model, type, dataset) %>%
  tally() %>%
  rowwise() %>%
  mutate(percent = n/100)

b02 <- stats_df %>%
  filter(type == "b0") %>%
  filter(dataset == "2") %>%
  group_by(prediction, model, type, dataset) %>%
  tally() %>%
  rowwise() %>%
  mutate(percent = n/100)

stats_df %>%
  filter(type == "b1") %>%
  #filter(dataset == "2") %>%
  group_by(prediction, model, type, dataset) %>%
  tally() %>%
  rowwise() %>%
  mutate(percent = n/100)

stats_df %>%
  filter(type == "b0") %>%
  mutate(prediction = factor(prediction, 
                             levels = c("incorrect", "correct"))) %>%
  ggplot(aes(x = mod, y = `50%`, 
             color = dataset)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_jitterdodge()) 

stats_df2 %>%
  filter(type == "b1") %>%
  mutate(prediction = factor(prediction, levels = c("incorrect", "correct")),
         var = factor(var, levels = c("low", "med", "high")),
         interval = factor(interval, levels = c("short", "long"))) %>%
  ggplot(aes(x = var, y = `50%`, 
             color = mod)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_jitterdodge()) +
  geom_hline(yintercept = b1) +
  facet_wrap(~interval)

ci_df <- stats_df2 %>%
  rowwise() %>%
  mutate(CI_width = `97.5%`-`2.5%`) %>%
  ungroup() 

ci_df %>%
  filter(type == "b1") %>%
  mutate(prediction = factor(prediction, levels = c("incorrect", "correct")),
         var = factor(var, levels = c("low", "med", "high")),
         interval = factor(interval, levels = c("short", "long"))) %>%
  ggplot(aes(x = var, y = CI_width, 
             fill = mod)) +
  geom_boxplot() +
  facet_wrap(~interval)

ci_df %>%
  filter(type == "b0") %>%
  mutate(prediction = factor(prediction, levels = c("incorrect", "correct")),
         var = factor(var, levels = c("low", "med", "high")),
         interval = factor(interval, levels = c("short", "long"))) %>%
  ggplot(aes(x = var, y = CI_width, 
             fill = mod)) +
  geom_boxplot() +
  facet_wrap(~interval)


stats_sum <- stats_df2 %>%
  group_by(type, model, interval, mod, var, prediction) %>%
  tally() %>%
  pivot_wider(names_from = "prediction",
              values_from = 'n') %>%
  rowwise() %>%
  mutate(prop_correct = correct/(correct + incorrect)) %>%
  ungroup() %>%
  mutate(prop_correct = case_when(correct == 100 ~ 1,
                                  TRUE ~ prop_correct)) 


stats_sum %>%
  mutate(var = factor(var, levels = c("low", "med", "high")),
         interval = factor(interval, levels = c("short", "long"))) %>%
ggplot(aes(x = var, y = prop_correct, fill = mod)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(type ~ interval)

