# Running the nest survival model
# Ana Miller-ter Kuile
# November 4, 2021

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

data <- readRDS(here("data_outputs", 
                     "01_whwonests",
                     '03_JAGS_input_data',
                     "mod3_JAGS_input_data.RDS"))

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
              "model3",
              "jags",
              "model3.R")

Sys.time() #~6 minutes
mod <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 4000,
                            DIC = TRUE)
Sys.time()

mcmcplot(mod$samples)
gelman.diag(mod$samples, multivariate = F)
# Check convergence -------------------------------------------------------

Rhat <- mod$Rhat
#both have converged
rhat_graph_fun2(Rhat)
rhat_graph_fun(Rhat)

gelman.diag(mod$samples, multivariate = F)

# Update for goodness of fit ----------------------------------------------

parms <- c("yrep_1", "yrep_2", 'resid_1', "resid_2", 'p1', "p2")

mod.update <- update(mod,
                     parameters.to.save = parms,
                     n.iter = 350,
                     codaOnly = TRUE)

mod_GOF <- mod.update

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- as.data.frame(data$y) %>%
  rename("y" = "data$y") %>%
  mutate(Nest_ID = 1:n()) 

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yrep_1 <- as.data.frame(mod_GOF$sims.list$yrep_1) %>%
  mutate(sample = 1:n()) %>%
  pivot_longer(1:(last_col()-1),
               values_to = "yrep_1",
               names_to = "Nest_ID") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, length(Nest_ID)))

yrep_2 <- as.data.frame(mod_GOF$sims.list$yrep_2) %>%
  mutate(sample = 1:n()) %>%
  pivot_longer(1:(last_col()-1),
               values_to = "yrep_2",
               names_to = "Nest_ID") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, length(Nest_ID)))

yrep <- yrep_1 %>%
  full_join(yrep_2, by = c("Nest_ID", "sample")) %>%
  mutate(yrep = case_when(!is.na(yrep_1) ~ yrep_1,
                                !is.na(yrep_2) ~ yrep_2,
                                TRUE ~ NA_real_)) %>%
  dplyr::select(-yrep_1, -yrep_2) %>%
  mutate(Nest_ID = as.numeric(Nest_ID))


# Confusion matrix dataset prep -------------------------------------------

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
                           'confusion_mod3.RDS'))

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


# Graph observed versus simulated -----------------------------------------
# 
# #posterior predictive check graphical observation
# (m3_pp <- ggplot() +
#    #graph the simulated data
#    geom_density(data = yrep, aes(x = Fate_class, group = Iteration, fill = type), 
#                 alpha = 0.2) +
#    geom_density(data = y, aes(x = Fate_class, fill = type), alpha = 0.5) + 
#    facet_wrap(~type))
# 
# #look like it varies by iteration how well it does
# 
# # Predicted accuracy ------------------------------------------------------
# 
# mu_p1 <- as.data.frame(mod_GOF$mean$yrep_1) %>%
#   rename("P" = 'mod_GOF$mean$yrep_1')
# 
# mu_p2 <- as.data.frame(mod_GOF$mean$yrep_2) %>%
#   rename("P" = 'mod_GOF$mean$yrep_2') %>%
#   filter(!is.na(P))
# 
# times <- as.data.frame(data$n.t) %>%
#   rename("intervals" = "data$n.t")
# 
# mu_p <- mu_p1 %>%
#   bind_rows(mu_p2)
# 
# y_acc <- y %>%
#   bind_cols(mu_p, times) %>%
#   mutate(Fate_class = as.factor(Fate_class)) %>%
#   mutate(type = case_when((Fate_class == 0 & P >= .5) ~ "mis-0",
#                           (Fate_class == 0 & P < .5) ~ "match-0",
#                           (Fate_class == 1 & P >= 0.5) ~ "match-1",
#                           (Fate_class == 1 & P < 0.5) ~ "mis-1"))
# 
# y_acc %>%
#   group_by(type) %>%
#   tally()
# 
# #accuary
# #0s:
# 80/(80+12) #87%
# #1s:
# 227/(227+11) #98%
# 
# (mod3_acc_plot <- ggplot(y_acc, aes(x = Fate_class, y = P)) +
#     geom_hline(yintercept = 0.5, linetype = 2) +
#     geom_boxplot() +
#     labs(x = "Observed fate",
#          y = "Predicted survival probability")  +
#     annotate(geom = "text", 
#              x = 0.75, y = 0.45,
#              label = "0%") +
#     annotate(geom = "text", 
#              x = 2.25, y = 0.55,
#              label = "95%"))
# 
# y_acc %>%
#   group_by(Fate_class) %>%
#   summarise(meanp = mean(P),
#             sdp = sd(P),
#             total = n(),
#             sep = sdp/sqrt(total))
# AUC ---------------------------------------------------------------------
resp <- as.vector(data$y)

iteration.num <- length(mod_GOF$sims.list$p1[,1])

pred1 <- as.data.frame(t(mod_GOF$sims.list$p1))
pred2 <- as.data.frame(t(mod_GOF$sims.list$p2))
pred2 <- pred2[29:nrow(pred2),]

pred <- rbind(pred1, pred2)

mod3_AUC <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod3_AUC[i] <- AUC_JAGS3(pred, 
                           iteration.num = i, 
                           resp = resp)
}

mean <- as.data.frame(mod3_AUC) %>%
  summarise(mean = mean(mod3_AUC)) %>%
  as_vector()

(mod3_AUC_plot <- as.data.frame(mod3_AUC) %>%
    ggplot() +
    geom_histogram(aes(x = mod3_AUC)) +
    geom_vline(xintercept = mean, linetype = 2) +
    labs(title = "Custom probability, logit link"))

saveRDS(mod3_AUC, 
        here('data_outputs',
             '01_whwonests',
             '05_for_plotting',
             'mod3AUC.RDS'))


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
             'Model3_posterior_samples.RDS'))
summary <- summary(mod$samples)

saveRDS(summary, 
        here('data_outputs',
             '01_whwonests',
             '04_posterior_summaries',
             'Model3_posterior_summary.RDS'))



