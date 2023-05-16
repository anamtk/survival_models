# Model pvalue results
# Ana Miller-ter Kuile
# May 15, 2023

# this script generates posterior plots for all model p-values


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


# Load posteriors ---------------------------------------------------------

m1 <- read.csv(here("data_outputs",
                    "02_kelp",
                    "04_posterior_summaries",
                    "Model1_pvalues.csv"))

m2 <- read.csv(here("data_outputs",
                    "02_kelp",
                    "04_posterior_summaries",
                    "Model2_pvalues.csv"))

m3 <- read.csv(here("data_outputs",
                    "02_kelp",
                    "04_posterior_summaries",
                    "Model3_pvalues.csv"))


# Merge data --------------------------------------------------------------

post <- bind_rows(m1, m2, m3)


# Filter out where all are NA ---------------------------------------------

post2 <- post %>%
  group_by(zvalue) %>%
  filter(any(p_cat == "s")) %>%
  mutate(group = case_when(zvalue %in% c("SST", "WavePower",
                                         "StipeNum") ~ "Interval",
                           TRUE ~ "Fixed")) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(mean_dir = Mean-0.5)

labs <- c("- \n p<0.05", "ns", "+ \n p<0.05")

post2 %>%
ggplot() +  
  geom_rect(aes(xmin = -Inf, 
                xmax = Inf, 
                ymin = -.45, 
                ymax = .45),
            fill = "lightgrey",
            alpha = 0.05) +
  geom_hline(yintercept = .45, linetype = 2, color = "lightgrey") +
  geom_hline(yintercept = -0.45, linetype = 2, color = "lightgrey") +
  # geom_bar(aes(x = zvalue, y = mean_dir,
  #                fill = Model),
  #            position = "dodge",
  #          stat = "identity",
  #          color = "gray47")+  
  geom_pointrange(aes(x = zvalue, y = mean_dir,
                      ymin= 0, ymax = mean_dir,
                 color = Model,shape = Model),
              size = 0.75,
             position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(labels = labs, breaks = c(-0.5, 0, 0.5)) +
  scale_color_manual(values = c('#8dd3c7','#bebada','#fdb462')) +
  coord_flip()  +
  labs(x = "Covariate", y = "Effect direction & Bayesian p-value") +
  facet_grid(group~., scales= "free", space = "free_y")

ggsave(filename = here("pictures",
                       "kelp_pvalues.pdf"),
       height = 6, 
       width = 7,
       units = "in")

