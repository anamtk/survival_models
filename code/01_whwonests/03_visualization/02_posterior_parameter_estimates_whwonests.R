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

theme_set(theme_bw())
# Load posterior summaries  -----------------------------------------------

mod1 <- readRDS(here('data_outputs',
                     '01_whwonests',
                     '04_posterior_summaries',
                     'Model1_posterior_summary.RDS'))

mod1daily <- readRDS(here('data_outputs',
                          '01_whwonests',
                          '04_posterior_summaries',
                          'Model1_Daily_posterior_summary.RDS'))

mod2 <- readRDS(here('data_outputs',
                     '01_whwonests',
                     '04_posterior_summaries',
                     'Model2_posterior_summary.RDS'))

mod2daily <- readRDS(here('data_outputs',
                          '01_whwonests',
                          '04_posterior_summaries',
                          'Model2daily_posterior_summary.RDS'))
mod2dailysubset <- readRDS(here('data_outputs',
                                '01_whwonests',
                                '04_posterior_summaries',
                                'Model2daily_subset_posterior_summary.RDS'))
mod3 <- readRDS(here('data_outputs',
                     '01_whwonests',
                     '04_posterior_summaries',
                     'Model3_posterior_summary.RDS'))

# Prep and combine --------------------------------------------------------

cov_effects_prep_fun <- function(model){
  
  modstats <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = 'parm') %>%
    filter(str_detect(parm, "b")) %>%
    filter(!str_detect(parm, "forest|year")) %>%
    mutate(Covariate = case_when(parm == "b0" ~ "Intercept",
                                 parm == 'b[1]' ~ "Nest height",
                                 parm == 'b[2]' ~ "Nest orientation",
                                 parm == 'b[3]' ~ "Initiation date",
                                 parm == 'b[4]' ~ "Small tree density",
                                 parm == 'b[5]' ~ "Large tree density",
                                 parm == 'b[6]' ~ "Maximum temperature",
                                 parm == 'b[7]' ~ "Maximum temperature^2",
                                 parm == 'b[8]' ~ "Precipitation",
                                 parm == 'b[9]' ~ "Precipitation^2",
                                 TRUE ~ NA_character_))
  
  return(modstats)
}


mod1stats <- cov_effects_prep_fun(model = mod1) %>%
  mutate(model = "Total")
mod1dailystats <- cov_effects_prep_fun(model = mod1daily) %>%
  mutate(model = "Total_Daily")
mod2stats <- cov_effects_prep_fun(model = mod2) %>%
  mutate(model = "Interval")
mod2dailystats <- cov_effects_prep_fun(model = mod2daily) %>%
  mutate(model = "Interval_Daily")
mod2dailysubsetstats <- cov_effects_prep_fun(model = mod2dailysubset) %>%
  mutate(model = "Interval_Daily_Subset")
mod3stats <- cov_effects_prep_fun(model = mod3) %>%
  mutate(model = "Custom")

all_stats <- bind_rows(mod1stats, mod1dailystats,
                       mod2stats, mod2dailystats,
                       mod2dailysubsetstats, mod3stats)


# Plot --------------------------------------------------------------------

all_stats <- all_stats %>%
  mutate(Covariate= factor(Covariate, levels = c("Maximum temperature^2",
                                                 "Maximum temperature",
                                                 "Precipitation^2",
                                                 "Precipitation",
                                                 "Small tree density",
                                                 "Large tree density",
                                                 "Nest height",          
                                                 "Nest orientation",
                                                 "Initiation date",
                                                 "Intercept"))) 

ggplot(all_stats, aes(x = Covariate, y = `50%`, color = model)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = `2.5%`,
                      ymax = `97.5%`),
                  position = position_dodge(0.9)) +
  coord_flip() +
  scale_color_manual(values = c('#9ebcda',
                                '#8c96c6',
                                '#8c6bb1',
                                '#88419d',
                                '#810f7c',
                                '#4d004b')) 

ggsave(here('pictures',
            'R',
            '01_whwonests',
            'nest_covariate_effects.jpg'),
       height = 4.5,
       width = 6,
       units = "in",
       dpi = 300)


# Custom vs interval ------------------------------------------------------


all_stats2 <- all_stats %>%
  filter(model %in% c("Total_Daily",
                      "Custom","Interval_Daily"))

covariate_effects_daily <- all_stats2 %>%
  #filter(Covariate != "Precipitation^2") %>%
ggplot(aes(x = `50%`, y = Covariate, color = model)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = `2.5%`,
                      xmax = `97.5%`),
                  position = position_dodge(0.9)) +
  #scale_x_sqrt() +
  scale_color_manual(values = c('#9ebcda',
                                '#8c96c6',
                                '#8c6bb1',
                                '#88419d',
                                '#810f7c')) 

