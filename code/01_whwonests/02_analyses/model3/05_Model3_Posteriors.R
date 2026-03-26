# Model covariate results
# Ana Miller-ter Kuile
# April 10, 2023

# this script generates posterior distributions
# for the survival model with interval-specific covariates
#but total end survival data


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

model3_sum <- readRDS(here('data_outputs',
                           '01_whwonests',
                           '04_posterior_summaries',
                           'Model3_posterior_summary.RDS'))


# Posterior median and CI of all parameters -------------------------------

parms <- c("b[1]", "b[2]",
           "b[3]", "b[4]","b[5]",            
           "b[6]", 'b[7]', 
           'b[8]', 'b[9]')

# "U", "B", "H", "HB"
#"PIPO", "Abies", "POTR5", "JUOC", "PSME"

mod3_est <- as.data.frame(model3_sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter %in% parms) %>%
  mutate(parameter = case_when(parameter == "b[1]" ~ "NestHt",
                               parameter == "b[2]" ~ "NestOrientation",
                               parameter == "b[3]" ~ "InitDay",
                               parameter == "b[4]" ~ "LgTreeDens",
                               parameter == "b[5]" ~ "SmTreeDens",
                               parameter == "b[6]" ~ "Tmax",
                               parameter == "b[7]" ~ "Tmax^2",
                               parameter == "b[8]" ~ "PPT",
                               parameter == "b[9]" ~ "PPT^2",
                               TRUE ~ parameter)) %>%
  mutate(Model = "Model3_CustomProb")

write.csv(mod3_est, here("data_outputs",
                         "01_whwonests",
                         "04_posterior_summaries",
                         "Model3_posteriors.csv"))



# Pvalues -----------------------------------------------------------------

p_values3 <- as.data.frame(model3_sum$statistics) %>%
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
  mutate(zvalue = case_when(zvalue == "z[1]" ~ "NestHt",
                            zvalue == "z[2]" ~ "NestOrientation",
                            zvalue == "z[3]" ~ "InitDay",
                            zvalue == "z[4]" ~ "LgTreeDens",
                            zvalue == "z[5]" ~ "SmTreeDens",
                            zvalue == "z[6]" ~ "Tmax",
                            zvalue == "z[7]" ~ "Tmax^2",
                            zvalue == "z[8]" ~ "PPT",
                            zvalue == "z[9]" ~ "PPT^2",
                            TRUE ~ zvalue)) %>%
  mutate(Model = "Model3_CustomProb")

p_values3 %>%
  filter(!is.na(p_dir)) %>%
  ggplot(aes(x = reorder(zvalue, p), y = p_dir)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  coord_flip()

write.csv(p_values3, here("data_outputs",
                          "01_whwonests",
                          "04_posterior_summaries",
                          "Model3_pvalues.csv"))
