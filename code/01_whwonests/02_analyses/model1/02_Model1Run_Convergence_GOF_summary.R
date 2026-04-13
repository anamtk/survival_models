# Get inits for survival Model 1
# Ana Miller-ter Kuile
# March 31, 2023

#this script runs an initial model for the full-survey interval 
#survival model and gets initials for monsoon runs

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  "R2jags", #jags wrapper
                  "coda",
                  "mcmcplots",
                  'pROC') #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#get model diagnostic plotting functions
source(here("code",
            "00_functions",
            "plot_functions.R"))

source(here("code",
            "00_functions",
            "GOF_functions.R"))

# Load Data ---------------------------------------------------------------
# 
# #load the formatted data for the JAGS model
data <- readRDS(here("data_outputs", 
                     "01_whwonests",
                     '03_JAGS_input_data',
                     "mod1_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------


params <- c(
            #Random covariate betas
            'b0.forest',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.forest',
            'sig.year',
            #covariates
            'b'
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "01_whwonests",
              "02_analyses",
              "model1",
              "jags",
              "model1.R")

mod <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 4000,
                            DIC = TRUE)

mcmcplot(mod$samples)

gelman.diag(mod$samples, multivariate = F)


# Graph Rhat --------------------------------------------------------------

Rhat <- mod$Rhat

#both have converged
rhat_graph_fun2(Rhat)
rhat_graph_fun(Rhat)

# Update for goodness of fit ----------------------------------------------

parms <- c("yrep", 'resid', 'p')

mod_GOF <- update(mod,
                     parameters.to.save = parms,
                     n.iter = 350,
                     codaOnly = TRUE)

saveRDS(mod_GOF, here("data_outputs",
                      '01_whwonests',
                      '04_posterior_summaries',
                      'mod1GOFsamples.RDS'))

mod_GOF <- readRDS(here("data_outputs",
                        '01_whwonests',
                        '04_posterior_summaries',
                        'mod1GOFsamples.RDS'))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  rename("y" = "data$y") %>%
  mutate(Nest_ID = 1:n()) 

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps <- mod_GOF$sims.list$yrep

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep<- reshape2::melt(yreps) %>%
  rename("sample" = "Var1",
         "Nest_ID" = "Var2",
         "yrep" = "value")

# prediction p ------------------------------------------------------------
p <- as.data.frame(mod_GOF$sims.list$p) %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-sample,
               names_to = "Nest_ID",
               values_to = "p_pos") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, nchar(Nest_ID))) %>%
  mutate(Nest_ID = as.numeric(Nest_ID),
         sample = as.numeric(sample)) %>%
  rowwise() %>%
  mutate(p_neg = 1-p_pos) %>%
  ungroup()


# Combine y, yrep, p and save ---------------------------------------------

predict_df <- yrep %>%
  left_join(p, by = c("Nest_ID", "sample")) %>%
  left_join(y, by = "Nest_ID")
  
saveRDS(predict_df, 
        here('data_outputs',
             '01_whwonests',
             '04_posterior_summaries',
             'mod1_yyrepp_df.RDS'))
# AUC ---------------------------------------------------------------------

resp <- as.vector(y$y)

iteration.num <- length(mod_GOF$sims.list$p[,1])

mod1_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod1_AUC[i] <- AUC_JAGS(mod_GOF, 
                          iteration.num = i, 
                          resp = resp)
}

mean <- as.data.frame(mod1_AUC) %>%
  summarise(mean = mean(mod1_AUC)) %>%
  as_vector()

(mod1_AUC_plot <- as.data.frame(mod1_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod1_AUC)) +
    geom_vline(xintercept = mean, linetype = 2) +
    labs(title = "Total survey exposure") )

saveRDS(mod1_AUC, 
        here('data_outputs',
             '01_whwonests',
             '05_for_plotting',
             'mod1AUC.RDS'))


# Confusion matrix --------------------------------------------------------

conf_df <- yrep %>%
  left_join(y, by = "Nest_ID")

conf_summary <- conf_df %>%
  mutate(category = case_when((yrep == 1 & y == 1) ~ "True Positive",
                              (yrep == 0 & y == 0) ~ "True Negative",
                              (yrep == 0 & y == 1) ~ "False Negative",
                              (yrep == 1 & y == 0) ~ "False Positive",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(category)) %>%
  group_by(sample, category) %>%
  tally() %>%
  ungroup() %>%
  mutate(Actually = case_when(category %in% c("True Positive",
                                              "False Negative") ~ "Actually Positive",
                              category %in% c("False Positive",
                                              "True Negative") ~ "Actually Negative",
                              TRUE ~ NA_character_),
         Predicted = case_when(category %in% c("True Positive",
                                               "False Positive") ~ "Predicted Positive",
                               category %in% c("False Negative",
                                               "True Negative") ~ "Predicted Negative"),
         Classification = case_when(str_detect(category, "True") ~ "Correct",
                                    str_detect(category, "False") ~ "Incorrect"))

ggplot(conf_summary, aes(x = 1, y = n, color = Classification)) +
  geom_boxplot() +
  facet_grid(Predicted ~ Actually, 
             scales = "free_x") +
  scale_y_log10() +
  scale_color_manual(values = c('#4d9221', 
                                '#c51b7d')) +
  theme(strip.background = element_blank(),
        axis.text.x  =element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  labs(y= "Count")

saveRDS(conf_summary, here("data_outputs", 
                           '01_whwonests',
                           '05_for_plotting',
                           'confusion_mod1.RDS'))

conf_sum <- conf_summary %>%
  ungroup() %>%
  group_by(category) %>%
  summarise(count = mean(n)) %>%
  pivot_wider(names_from = 'category',
              values_from = 'count')

conf_sum %>%
  reframe(Accuracy = (`True Positive` + `True Negative`)/
            (`True Positive` + `True Negative` + `False Positive` + `False Negative`),
          Balanced_Accuracy = ((`True Positive`/(`True Positive` + `False Negative`)) +
            (`True Negative`/(`True Negative` + `False Positive`)))/2)

# Save summaries ----------------------------------------------------------

# Output summary stats and coda samples -----------------------------------

samples <- mod$samples

#pull all three chains out of the mcmc object
a <- as.matrix(samples[[1]])
b <- as.matrix(samples[[2]])
c <- as.matrix(samples[[3]])

#make this a dataframe bound by rows
samples_all <- as.data.frame(rbind(a,b,c))

#subset just a few for now because this is a LOT - can change the n 
# later if we want - I would probably aim for ~1000+ for any final
# analyses
samples_all2 <- slice_sample(samples_all, n = 4000)

saveRDS(samples_all2, 
        here('data_outputs',
             '01_whwonests',
             '04_posterior_summaries',
             'Model1_posterior_samples.RDS'))
        
summary <- summary(mod$samples)

saveRDS(summary, 
        here('data_outputs',
             '01_whwonests',
             '04_posterior_summaries',
             'Model1_posterior_summary.RDS'))



