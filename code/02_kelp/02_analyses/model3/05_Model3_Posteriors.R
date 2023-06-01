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

model3_sum <- readRDS(here("monsoon",
                           "02_kelp",
                           "model3",
                           "outputs",
                           "model3_posterior_summary.RDS"))


# Posterior median and CI of all parameters -------------------------------

parms <- c("b[1]", "b[2]",
           "b[3]", "b[4]","b[5]", 
           'b[6]', 'b[7]', 'b[8]',
           'b9[2]', 'b9[3]',
           'b9[4]', 'b9[5]')

# "B", "BO",
# "C", "S", "SS"

mod3_est <- as.data.frame(model3_sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter %in% parms) %>%
  mutate(parameter = case_when(parameter == "b[1]" ~ "Sea surface temp",
                               parameter == "b[2]" ~ "Wave power",
                               parameter == "b[3]" ~ "Stipe number",
                               parameter == "b[4]" ~ "Holdfast diameter",
                               parameter == "b[5]" ~ "Plant depth",
                               parameter == "b[6]" ~ "Wave power*Depth",
                               parameter == "b[7]" ~ "Wave power*Stipe number",
                               parameter == "b[8]" ~ "Wave power*Holdfast diameter",
                               parameter == 'b9[2]' ~ "Substrate = boulder",
                               parameter == "b9[3]" ~ "Substrate = cobble",
                               parameter == "b9[4]" ~ "Substrate = sand",
                               parameter == "b9[5]" ~ "Substrate = shallow sand",
                               TRUE ~ parameter)) %>%
  mutate(Model = "Model3_CustomProb")

write.csv(mod3_est, here("data_outputs",
                         "02_kelp",
                         "04_posterior_summaries",
                         "Model3_posteriors.csv"))


# Pvalues -----------------------------------------------------------------

# "B", "BO",
# "C", "S", "SS"

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
  filter(!zvalue %in% c("z.b9[1]")) %>%
  mutate(zvalue = case_when(zvalue == 'z.b9[2]' ~ 
                              'Substrate:Boulder',
                            zvalue == "z.b9[3]" ~
                              "Substrate:Cobble",
                            zvalue == "z.b9[4]" ~
                              "Substrate:Sand",
                            zvalue == "z.b9[5]" ~
                              "Substrate:ShallowSand",
                            zvalue == "z[1]" ~ "SST",
                            zvalue == "z[2]" ~ "WavePower",
                            zvalue == "z[3]" ~ "StipeNum",
                            zvalue == "z[4]" ~ "HoldfastDiam",
                            zvalue == "z[5]" ~ "Depth",
                            zvalue == "z[6]" ~ "Wave power*Depth",
                            zvalue == "z[7]" ~ "Wave power*Stipe number",
                            zvalue == "z[8]" ~ "Wave power*Holdfast diameter",
                            TRUE ~ zvalue)) %>%
  mutate(Model = "Model3_CustomProb")

p_values3 %>%
  filter(!is.na(p_dir)) %>%
  ggplot(aes(x = reorder(zvalue, p), y = p_dir)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  coord_flip()

write.csv(p_values3, here("data_outputs",
                          "02_kelp",
                          "04_posterior_summaries",
                          "Model3_pvalues.csv"))



#