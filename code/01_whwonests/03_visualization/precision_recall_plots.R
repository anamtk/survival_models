
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "lubridate", "glmmTMB",
                  'yardstick')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load gof datasets -------------------------------------------------------

mod1 <- readRDS(here('data_outputs',
                          '01_whwonests',
                          '04_posterior_summaries',
                          'mod1_yyrepp_df.RDS')) %>%
  mutate(model = "total")

mod1daily <- readRDS(here('data_outputs',
                     '01_whwonests',
                     '04_posterior_summaries',
                     'mod1daily_yyrepp_df.RDS'))%>%
  mutate(model = "totaldaily")

mod2 <- readRDS(here('data_outputs',
                          '01_whwonests',
                          '04_posterior_summaries',
                          'mod2_yyrepp_df.RDS'))%>%
  mutate(model = "interval")

mod2daily <- readRDS(here('data_outputs',
                          '01_whwonests',
                          '04_posterior_summaries',
                          'mod2daily_yyrepp_df.RDS'))%>%
  mutate(model = "intervaldaily")

mod3 <- readRDS(here('data_outputs',
                          '01_whwonests',
                          '04_posterior_summaries',
                          'mod3_yyrepp_df.RDS'))%>%
  mutate(model = "custom")


# Combine all of them -----------------------------------------------------

all_predictions <- mod1 %>%
  bind_rows(mod1daily, mod2, mod2daily, mod3) %>%
  dplyr::select(Nest_ID, y, yrep, p_pos, p_neg, model, sample)

# Compute precision and recall --------------------------------------------

predict_summary <- all_predictions %>%
  mutate(category = case_when((yrep == 1 & y == 1) ~ "True Positive",
                              (yrep == 0 & y == 0) ~ "True Negative",
                              (yrep == 0 & y == 1) ~ "False Negative",
                              (yrep == 1 & y == 0) ~ "False Positive",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(category)) %>%
  group_by(sample, model, category) %>%
  tally() %>%
  ungroup() 

conf_sum <- predict_summary %>%
  ungroup() %>%
  group_by(category, model) %>%
  summarise(count = mean(n)) %>%
  pivot_wider(names_from = 'category',
              values_from = 'count')

conf_sum %>%
  group_by(model) %>%
  reframe(Accuracy = (`True Positive` + `True Negative`)/
            (`True Positive` + `True Negative` + `False Positive` + `False Negative`),
          Balanced_Accuracy = ((`True Positive`/(`True Positive` + `False Negative`)) +
                                 (`True Negative`/(`True Negative` + `False Positive`)))/2,
          Neg_Recall = `True Negative`/(`True Negative` + `False Positive`),
          Neg_Precision = `True Negative`/(`True Negative` + `False Negative`),
          Neg_F1 = (2*(Neg_Recall*Neg_Precision)/(Neg_Precision + Neg_Recall))
  )


# Get PR-AUC --------------------------------------------------------------

all_predictions %>%
  mutate(y = factor(y, levels = c("0", "1"))) %>%
  group_by(model, sample) %>%
  pr_auc(truth = y, p_neg) %>%
  ggplot(aes(x = .estimate)) +
  geom_histogram()+
  facet_wrap(~model, scales= 'free_y')

subsample <- sample(mod1$sample, size = 10,
                    replace = F)

cutoffs <- as.data.frame(bind_cols(start = seq(0,0.98, by = 0.02),
                                   end = seq(0.02, 1, by = 0.02))) %>%
  mutate(group = 1:n()) 

pr_curves <- all_predictions %>%
  mutate(y = factor(y, levels = c("0", "1"))) %>%
  group_by(model, sample) %>%
  pr_curve(truth = y, p_neg) %>%
  filter(.threshold < Inf & .threshold > -Inf) %>%
  ungroup() %>%
  left_join(cutoffs, 
            by = join_by(between(.threshold, start, end))) %>%
  group_by(model, group, start) %>%
  reframe(mean_recall = mean(recall),
          mean_precision = mean(precision),
          sd_recall = sd(recall),
          sd_precision = sd(precision)) %>%
  ungroup()

pr_curves_daily <- pr_curves %>%
  filter(model %in% c('custom',
                      'intervaldaily',
                      'totaldaily')) %>%
ggplot(aes(x = mean_recall, y = mean_precision)) +
  geom_point(aes(color = model)) +
  geom_linerange(aes(xmin = mean_recall -sd_recall,
                    xmax = mean_recall + sd_recall,
                    color = model)) +
  geom_linerange(aes(ymin = mean_precision - sd_precision,
                    ymax = mean_precision + sd_precision,
                    color= model)) +
  scale_color_manual(values = c('#9ebcda',
                                '#8c96c6',
                                '#8c6bb1',
                                '#88419d',
                                '#810f7c')) +
  labs(x = "Recall \n (how many actual negatives did model identify)", 
       y = "Precision \n (how many predicted negatives \n are actually negative)")


