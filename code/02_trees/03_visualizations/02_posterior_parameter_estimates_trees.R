#covariate effect plots for the
#nest dataset

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


# Load posterior summaries  -----------------------------------------------

mod1 <- readRDS(here('data_outputs',
                     '02_trees',
                     '04_posterior_summaries',
                     'model1',
                     'summary_model1.RDS'))

mod2 <- readRDS(here('data_outputs',
                     '02_trees',
                     '04_posterior_summaries',
                     'model2',
                     'summary_model2.RDS'))

mod3 <- readRDS(here('data_outputs',
                    '02_trees',
                    '04_posterior_summaries',
                    'model3',
                    'summary_model3.RDS'))

# Prep and combine --------------------------------------------------------

cov_effects_prep_fun <- function(model){
  
  modstats <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = 'parm') %>%
    filter(str_detect(parm, "b")) %>%
    filter(!str_detect(parm, "forest|year")) %>%
    mutate(Covariate = case_when(parm == "b0" ~ "Intercept",
                                 parm == 'b[2]' ~ "DBH",
                                 parm == 'b[3]' ~ "DBH^2",
                                 parm == 'b[4]' ~ "Basal area",
                                 parm == 'b[5]' ~ "Canopy cover",
                                 parm == 'b[6]' ~ "Maximum VPD",
                                 parm == 'b[7]' ~ "Minimum SWA",
                                 parm == 'b1[2]' ~ "Treated",
                                 TRUE ~ NA_character_))
  
  return(modstats)
}


mod1stats <- cov_effects_prep_fun(model = mod1) %>%
  mutate(model = "Total")
mod2stats <- cov_effects_prep_fun(model = mod2) %>%
  mutate(model = "Interval")
mod3stats <- cov_effects_prep_fun(model = mod3) %>%
  mutate(model = "Custom")

all_stats <- bind_rows(mod1stats, mod2stats, mod3stats) %>%
  filter(parm != "b1[1]")


# Plot --------------------------------------------------------------------

all_stats <- all_stats %>%
mutate(Covariate= factor(Covariate, levels = c("Maximum VPD",
                                               "Minimum SWA",
                                               "Canopy cover",
                                               "Basal area",
                                               "DBH^2",
                                               "DBH",
                                               "Treated",
                                               "Intercept")))

ggplot(all_stats, aes(x = Covariate, y = `50%`, color = model)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = `2.5%`,
                      ymax = `97.5%`),
                  position = position_dodge(0.9)) +
  coord_flip() +
  scale_color_manual(values = c('#8c96c6',
                                '#8c6bb1',
                                '#88419d')) 

ggsave(here('pictures',
            'R',
            '02_trees',
            'tree_covariate_effects.jpg'),
       height = 4,
       width = 4.5,
       units = "in",
       dpi = 300)


# Interval vs custom ------------------------------------------------------


all_stats2 <- all_stats %>%
  filter(model != "Total")

ggplot(all_stats2, aes(x = Covariate, y = `50%`, color = model)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = `2.5%`,
                      ymax = `97.5%`),
                  position = position_dodge(0.9)) +
  coord_flip() +
  scale_color_manual(values = c('#8c96c6',
                                '#8c6bb1')) 
