#AUC for the nest dataset

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "patchwork") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load AUC datasets -------------------------------------------------------

mod1AUC <- readRDS(here('data_outputs',
                        '01_whwonests',
                        '05_for_plotting',
                        'mod1AUC.RDS')) %>%
  as.data.frame() %>%
  rename("AUC" = ".") %>%
  mutate(model = "Total")

mod2AUC <- readRDS(here('data_outputs',
                        '01_whwonests',
                        '05_for_plotting',
                        'mod2AUC.RDS'))%>%
  as.data.frame() %>%
  rename("AUC" = ".") %>%
  mutate(model = "Interval")

mod3AUC <- readRDS(here('data_outputs',
                        '01_whwonests',
                        '05_for_plotting',
                        'mod3AUC.RDS'))%>%
  as.data.frame() %>%
  rename("AUC" = ".") %>%
  mutate(model = "Custom")


# Combine -----------------------------------------------------------------

allAUC <- bind_rows(mod1AUC, mod2AUC, mod3AUC)


# Get means ---------------------------------------------------------------

means <- allAUC %>%
  group_by(model) %>%
  summarise(mean = mean(AUC)) %>%
  ungroup()


# Plot histogram ----------------------------------------------------------

ggplot() +
  geom_histogram(data = allAUC, 
                 aes(x = AUC, 
                     fill = model), 
                 alpha = 0.8) +
  geom_vline(data = means, 
             aes(xintercept = mean,
                 color = model),
             linetype = 2) +
  scale_fill_manual(values = c('#8c96c6',
                               '#8c6bb1',
                               '#88419d')) +
  scale_color_manual(values = c('#8c96c6',
                               '#8c6bb1',
                               '#88419d')) +
  labs(y = "Number of posterior samples")


ggsave(here('pictures',
            'R',
            '01_whwonests',
            'nest_AUC.jpg'),
       height = 3,
       width = 4.5,
       units = "in",
       dpi = 300)
# Confusion matrices ------------------------------------------------------


# Load confusion datasets -------------------------------------------------

conf1 <- readRDS(here("data_outputs", 
                      '01_whwonests',
                      '05_for_plotting',
                      'confusion_mod1.RDS')) %>%
  dplyr::select(-sample) %>%
  mutate(model = 'Total')


conf2 <- readRDS(here("data_outputs", 
                      '01_whwonests',
                      '05_for_plotting',
                      'confusion_mod2.RDS'))%>%
  dplyr::select(-sample) %>%
  mutate(model = 'Interval')


conf3 <- readRDS(here("data_outputs", 
                      '01_whwonests',
                      '05_for_plotting',
                      'confusion_mod3.RDS'))%>%
  dplyr::select(-sample) %>%
  mutate(model = 'Custom')


# Combine confusion -------------------------------------------------------

conf_all <- bind_rows(conf1, conf2, conf3)

ggplot(conf_all, aes(x = model, y = n, color = Classification)) +
  geom_boxplot() +
  facet_grid(Predicted ~ Actually, 
             scales = "free_x") +
  scale_y_log10() +
  scale_color_manual(values = c('#4d9221', 
                                '#c51b7d')) +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  labs(y= "Count")

# Accuracy-Balanced Accuracy ----------------------------------------------

#adapt this by model:
conf_sum <- conf_all %>%
  ungroup() %>%
  group_by(category, model) %>%
  summarise(count = mean(n)) %>%
  pivot_wider(names_from = 'category',
              values_from = 'count')

conf_sum %>%
  ungroup() %>%
  group_by(model) %>%
  #The Neg_Precision is the fraction of True - elements divided 
  #by the total number of - predicted units
  #The Neg_Recall is the fraction of True - elements 
  #divided by the total number of - classified units
  reframe(Accuracy = (`True Positive` + `True Negative`)/
            (`True Positive` + `True Negative` + `False Positive` + `False Negative`),
          Balanced_Accuracy = ((`True Positive`/(`True Positive` + `False Negative`)) +
                                 (`True Negative`/(`True Negative` + `False Positive`)))/2,
          Neg_Recall = `True Negative`/(`True Negative` + `False Positive`),
          Neg_Precision = `True Negative`/(`True Negative` + `False Negative`),
          Neg_F1 = (2*(Neg_Recall*Neg_Precision)/(Neg_Precision + Neg_Recall))
          )

