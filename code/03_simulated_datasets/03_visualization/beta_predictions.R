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

b0 <- 2.8
b1 <- 3.5

# Low var -----------------------------------------------------------------


# Load datasets -----------------------------------------------------------

mod1low <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod1low_summary.RDS'))

mod2low <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod2low_summary.RDS'))

mod3low <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod3low_summary.RDS'))

# Prep data ---------------------------------------------------------------


lowvar <- post_fun(mod1 = mod1low,
         mod2 = mod2low,
         mod3 = mod3low) %>%
  mutate(model = case_when(model == "Total Exposure" ~ "Total",
                           TRUE ~ model))

(lowb0 <- ggplot(lowvar, aes(x = model, y = b0)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.2) +
  geom_hline(yintercept = b0, linetype = "dashed") +
  labs(title = "Low covariate variation",
       y = "Intercept")  +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylim(-9, 9))

(lowb1 <- ggplot(lowvar, aes(x = model, y = b1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.2) +
  geom_hline(yintercept = b1, linetype = "dashed") +
  labs(y = "Covariate effect")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  ylim(-25, 28))


# Med var -----------------------------------------------------------------

# Load datasets -----------------------------------------------------------

mod1med <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod1med_summary.RDS'))

mod2med <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod2med_summary.RDS'))

mod3med <- readRDS(here('data_outputs',
                        '03_simulated',
                        '04_posterior_summaries',
                        'mod3med_summary.RDS'))

# Prep data ---------------------------------------------------------------

medvar <- post_fun(mod1 = mod1med,
                   mod2 = mod2med,
                   mod3 = mod3med)%>%
  mutate(model = case_when(model == "Total Exposure" ~ "Total",
                           TRUE ~ model))

(medb0 <- ggplot(medvar, aes(x = model, y = b0)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.2) +
  geom_hline(yintercept = b0, linetype = "dashed") +
  labs(title = "Medium covariate variation",
       y = "Intercept")  +
  theme(axis.text = element_blank(),
        axis.title = element_blank()) +
  ylim(-9, 9))

(medb1 <- ggplot(medvar, aes(x = model, y = b1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.2) +
  geom_hline(yintercept = b1, linetype = "dashed") +
  labs(y = "Covariate effect")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(-25, 28))


# High var ----------------------------------------------------------------

# Load datasets -----------------------------------------------------------

mod1high <- readRDS(here('data_outputs',
                         '03_simulated',
                         '04_posterior_summaries',
                         'mod1high_summary.RDS'))

mod2high <- readRDS(here('data_outputs',
                         '03_simulated',
                         '04_posterior_summaries',
                         'mod2high_summary.RDS'))

mod3high <- readRDS(here('data_outputs',
                         '03_simulated',
                         '04_posterior_summaries',
                         'mod3high_summary.RDS'))

# Prep data ---------------------------------------------------------------

highvar <- post_fun(mod1 = mod1high,
                   mod2 = mod2high,
                   mod3 = mod3high)%>%
  mutate(model = case_when(model == "Total Exposure" ~ "Total",
                           TRUE ~ model))

(highb0 <- ggplot(highvar, aes(x = model, y = b0)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.2) +
  geom_hline(yintercept = b0, linetype = "dashed") +
  labs(title = "High covariate variation",
       y = "Intercept")  +
  theme(axis.text = element_blank(),
        axis.title = element_blank()) +
  ylim(-9, 9))

(highb1 <- ggplot(highvar, aes(x = model, y = b1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.2) +
  geom_hline(yintercept = b1, linetype = "dashed") +
  labs(y = "Covariate effect") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.title = element_blank()) +
  ylim(-25, 28) )



# Plot together -----------------------------------------------------------

p1 <- lowb0 + medb0 + highb0

p2 <- lowb1 + medb1 + highb1

(p3 <- (lowb0 + medb0 + highb0)/(lowb1 + medb1 + highb1))


# Remove total exposure and plot ------------------------------------------

# lowvar2 <- lowvar %>%
#   filter(model != "Total Exposure")
#   
# medvar2 <- medvar %>%
#   filter(model != "Total Exposure")
# 
# highvar2 <- highvar %>%
#   filter(model != "Total Exposure")
# 
# lowb02 <- ggplot(lowvar2, aes(x = model, y = b0)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_hline(yintercept = b0, linetype = "dashed") +
#   labs(title = "Low covariate variation")
# 
# lowb12 <- ggplot(lowvar2, aes(x = model, y = b1)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_hline(yintercept = b1, linetype = "dashed") +
#   labs(title  = "Low covariate variation")
# 
# medb02 <- ggplot(medvar2, aes(x = model, y = b0)) +
#   geom_boxplot() +
#   geom_hline(yintercept = b0, linetype = "dashed") +
#   labs(title = "Medium covariate variation")  
# 
# medb12 <- ggplot(medvar2, aes(x = model, y = b1)) +
#   geom_boxplot() +
#   geom_hline(yintercept = b1, linetype = "dashed") +
#   labs(title = "Medium covariate variation")
# 
# highb02 <- ggplot(highvar2, aes(x = model, y = b0)) +
#   geom_boxplot() +
#   geom_hline(yintercept = b0, linetype = "dashed") +
#   labs(title = "High covariate variation")
# 
# highb12 <- ggplot(highvar2, aes(x = model, y = b1)) +
#   geom_boxplot() +
#   geom_hline(yintercept = b1, linetype = "dashed") +
#   labs(title = "High covariate variation")
# 
# p12 <- lowb02 + medb02 + highb02
# 
# p22 <- lowb12 + medb12 + highb12
# 
# (p32 <- (lowb02 + medb02 + highb02)/(lowb12 + medb12 + highb12))


