# Model covariate results
# Ana Miller-ter Kuile
# April 10, 2023

# this script generates posterior distributions
# for the survival model with full interval


# Load packages -----------------------------------------------------------


# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "patchwork")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here::here("code", 
                  "00_functions",
                  "plot_functions.R"))

source(here::here("code",
                  "00_functions",
                  "tidy_functions.R"))

# Load data ---------------------------------------------------------------

model2_sum <- readRDS(here("monsoon",
                           "03_trees",
                           "model2",
                           "outputs",
                           "model2_posterior_summary.RDS"))


# Posterior median and CI of all parameters -------------------------------

parms <- c("b1[1]", "b1[2]", 'b[2]',
           "b[3]", "b[4]","b[5]",  
           'b[6]', 'b[7]')


mod2_est <- as.data.frame(model2_sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter %in% parms) %>%
  mutate(parameter = case_when(parameter == "b1[1]" ~ "Treatment: 1",
                               parameter == "b1[2]" ~ "Treatment: 2",
                               parameter == "b[2]" ~ "DBH",
                               parameter == "b[3]" ~ "DBH^2",
                               parameter == "b[4]" ~ "Basal area",
                               parameter == "b[5]" ~ "Canopy cover",
                               parameter == "b[6]" ~ "VPD",
                               parameter == "b[7]" ~ "SWA",
                               TRUE ~ parameter)) %>%
  mutate(Model = "Model2_IntervalData")
  
write.csv(mod2_est, here("data_outputs",
                         "03_trees",
               "04_posterior_summaries",
               "Model2_posteriors.csv"))


# Pvalues -----------------------------------------------------------------

# "B", "BO",
# "C", "S", "SS"

p_values2 <- as.data.frame(model2_sum$statistics) %>%
  dplyr::select(Mean) %>%
  rownames_to_column(var = "zvalue") %>%
  filter(str_detect(zvalue, "z")) %>%
  mutate(p = case_when(Mean >= 0.5 ~ (1-Mean), #I Think these are 1-tailed p-values but check with Kiona
                       Mean < 0.5 ~ (1 - (1-Mean)))) %>%
  mutate(p_cat = case_when(p <= 0.05 ~ "s",
                           TRUE ~ "ns")) %>%
  mutate(direction = case_when(Mean >= 0.5 ~ 'positive', 
                               Mean < 0.5 ~ "negative")) %>%
  mutate(p_dir = case_when((p_cat == "s" & direction == "positive") ~ p,
                           (p_cat == "s" & direction == "negative") ~ -p,
                           TRUE ~ NA_real_)) %>%
  filter(!zvalue %in% c("z.b6[1]")) %>%
  mutate(zvalue = case_when(zvalue == 'z.b6[2]' ~ 
                              'Substrate:Boulder',
                            zvalue == "z.b6[3]" ~
                              "Substrate:Cobble",
                            zvalue == "z.b6[4]" ~
                              "Substrate:Sand",
                            zvalue == "z.b6[5]" ~
                              "Substrate:ShallowSand",
                            zvalue == "z[1]" ~ "SST",
                            zvalue == "z[2]" ~ "WavePower",
                            zvalue == "z[3]" ~ "StipeNum",
                            zvalue == "z[4]" ~ "HoldfastDiam",
                            zvalue == "z[5]" ~ "Depth",
                            TRUE ~ zvalue)) %>%
  mutate(Model = "Model2_IntervalData")

p_values2 %>%
  filter(!is.na(p_dir)) %>%
  ggplot(aes(x = reorder(zvalue, p), y = p_dir)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  coord_flip()

write.csv(p_values2, here("data_outputs",
                          "02_kelp",
                          "04_posterior_summaries",
                          "Model2_pvalues.csv"))




