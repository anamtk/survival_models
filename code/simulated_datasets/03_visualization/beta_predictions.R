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


# Low var -----------------------------------------------------------------


# Load datasets -----------------------------------------------------------

mod1low <- readRDS(here("monsoon",
                     "simulated",
                     "lowvar",
                     "outputs",
                     'mod1_lowvar_sum.RDS'))

mod2low <- readRDS(here("monsoon",
                     "simulated",
                     "lowvar",
                     "outputs",
                     'mod2_lowvar_sum.RDS'))

mod3low <- readRDS(here("monsoon",
                     "simulated",
                     "lowvar",
                     "outputs",
                     'mod3_lowvar_sum.RDS'))

# Prep data ---------------------------------------------------------------

#b0 <- 0.05
#b1 <- 0.5 #this value makes sure that ps is always >0.5 in the range of x

lowvar <- post_fun(mod1 = mod1low,
         mod2 = mod2low,
         mod3 = mod3low)

lowb0 <- ggplot(lowvar, aes(x = model, y = b0)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #geom_hline(yintercept = 1.38, linetype = "dashed") +
  labs(title = "Low covariate variation") +
  ylim(-17,2)

lowb1 <- ggplot(lowvar, aes(x = model, y = b1)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #geom_hline(yintercept = 1.5, linetype = "dashed") +
  labs(title  = "Low covariate variation") +
  ylim(-0.1, 6)


# Med var -----------------------------------------------------------------

# Load datasets -----------------------------------------------------------

mod1med <- readRDS(here("monsoon",
                        "simulated",
                        "medvar",
                        "outputs",
                        'mod1_medvar_sum.RDS'))

mod2med <- readRDS(here("monsoon",
                        "simulated",
                        "medvar",
                        "outputs",
                        'mod2_medvar_sum.RDS'))

mod3med <- readRDS(here("monsoon",
                        "simulated",
                        "medvar",
                        "outputs",
                        'mod3_medvar_sum.RDS'))

# Prep data ---------------------------------------------------------------

#b0 <- 1.38 #gets mean value survival fairly high
#b1 <- 0.5 #this value makes sure that ps is always >0.5 in the range of x

medvar <- post_fun(mod1 = mod1med,
                   mod2 = mod2med,
                   mod3 = mod3med)

medb0 <- ggplot(medvar, aes(x = model, y = b0)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #geom_hline(yintercept = 1.38, linetype = "dashed") +
  labs(title = "Medium covariate variation") +
  ylim(-17,2)

medb1 <- ggplot(medvar, aes(x = model, y = b1)) +
  geom_boxplot() +
  #geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = "Medium covariate variation")+
  ylim(-0.1, 6)


# High var ----------------------------------------------------------------

# Load datasets -----------------------------------------------------------

mod1high <- readRDS(here("monsoon",
                        "simulated",
                        "highvar",
                        "outputs",
                        'mod1_highvar_sum.RDS'))

mod2high <- readRDS(here("monsoon",
                        "simulated",
                        "highvar",
                        "outputs",
                        'mod2_highvar_sum.RDS'))

mod3high <- readRDS(here("monsoon",
                        "simulated",
                        "highvar",
                        "outputs",
                        'mod3_highvar_sum.RDS'))

# Prep data ---------------------------------------------------------------

#b0 <- 1.38 #gets mean value survival fairly high
#b1 <- 0.5 #this value makes sure that ps is always >0.5 in the range of x

highvar <- post_fun(mod1 = mod1high,
                   mod2 = mod2high,
                   mod3 = mod3high)

highb0 <- ggplot(highvar, aes(x = model, y = b0)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #geom_hline(yintercept = 1.38, linetype = "dashed") +
  labs(title = "High covariate variation") +
  ylim(-17,2)

highb1 <- ggplot(highvar, aes(x = model, y = b1)) +
  geom_boxplot() +
  #geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = "High covariate variation")+
  ylim(-0.1, 6)



# Plot together -----------------------------------------------------------

lowb0 + medb0 + highb0

lowb1 + medb1 + highb1

(lowb0 + medb0 + highb0)/(lowb1 + medb1 + highb1)
