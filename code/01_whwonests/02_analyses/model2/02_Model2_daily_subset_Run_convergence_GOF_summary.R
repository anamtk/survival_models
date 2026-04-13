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
                     "mod2_daily_subset_JAGS_input_data.RDS"))

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
              "model2",
              "jags",
              "model2_interval_daily.R")

Sys.time() #~ minutes
mod <- jagsUI::jags(data = data,
                            inits = NULL,
                            model.file = model,
                            parameters.to.save = params,
                            parallel = TRUE,
                            n.chains = 3,
                            n.iter = 5000,
                    n.burnin = 1000,
                            DIC = TRUE)
Sys.time()
mcmcplot(mod$samples)

# Check convergence -------------------------------------------------------

Rhat <- mod$Rhat
#both have converged
rhat_graph_fun2(Rhat)
rhat_graph_fun(Rhat)

gelman.diag(mod$samples, multivariate = F)


# GOF ---------------------------------------------------------------------

parms <- c("y.repkeep", 'resid', 'p.intkeep', 'p.int', 'yrep')


mod_GOF <- update(mod,
                  parameters.to.save = parms,
                  n.iter = 350,
                  codaOnly = TRUE)

# AUC ---------------------------------------------------------------------

#resp <- as.vector(y$Fate_class)

# iteration.num <- length(mod_GOF$sims.list$p.intkeep[,1])
# 
# AUC_JAGS4(mod_GOF, 
#           iteration.num = 11, 
#           resp = resp)
# 
# mod2_AUC <- rep(NA, iteration.num)
# 
# for(i in 1:iteration.num){
#   mod2_AUC[i] <- AUC_JAGS4(mod_GOF, 
#                            iteration.num = i, 
#                            resp = resp)
# }
# 
# mean <- as.data.frame(mod2_AUC) %>%
#   summarise(mean = mean(mod2_AUC)) %>%
#   as_vector()
# 
# (mod2_AUC_plot <- as.data.frame(mod2_AUC) %>%
#     ggplot() +
#     geom_histogram(aes(x = mod2_AUC)) +
#     geom_vline(xintercept = mean, linetype = 2) +
#     labs(title = "Interval-level response \n (last survey data only), logit link"))


#THIS IS for all - would need to track pint though for all intervals
t <- mod_GOF$sims.list$p.int

layers <- dim(t)[[3]]

dfs <- lapply(1:layers,
              function(x){
                return(as.data.frame(t[,,x]))
              } )

dfs1 <- dfs %>%
  map(~mutate(., iteration = 1:n()))

full_df <- bind_rows(dfs1, .id = "interval") %>%
  pivot_longer(cols = 2:(last_col()-1),
               values_to = "p",
               names_to = "Nest_ID") %>%
  mutate(Nest_ID = str_sub(Nest_ID, 2, length(Nest_ID))) %>%
  unite(col = "ID_interval",
        c("Nest_ID", "interval"),
        sep = "_") %>%
  filter(!is.na(p))


resp <- as.data.frame(data$y) %>%
  mutate(Nest_ID = 1:n()) %>%
  pivot_longer(cols = 1:(last_col()-1),
               values_to = "resp",
               names_to = "interval") %>%
  filter(!is.na(resp)) %>%
  unite(col = "ID_interval",
        c("Nest_ID", "interval"),
        sep = "_")

iteration.num <- length(unique(full_df$iteration))

mod2_AUC2 <- rep(NA, iteration.num)

for(i in 1:iteration.num){
  mod2_AUC2[i] <- AUC_JAGS2(df = full_df,
                            iteration.num = i,
                            resp = resp$resp)
}

mean2 <- as.data.frame(mod2_AUC2) %>%
  summarise(mean = mean(mod2_AUC2)) %>%
  as_vector()

(mod2_AUC_plotall <- as.data.frame(mod2_AUC2) %>%
    ggplot() +
    geom_histogram(aes(x = mod2_AUC2)) +
    geom_vline(xintercept = mean2, linetype = 2) +
    labs(x = "AUC",
         y = "Number of samples"))

saveRDS(mod2_AUC2, 
        here('data_outputs',
             '01_whwonests',
             '05_for_plotting',
             'mod2dailyAUC.RDS'))

# Confusion matrix --------------------------------------------------------

#THIS IS for all - would need to track pint though for all intervals
yrep <- mod_GOF$sims.list$yrep

dimnames(yrep) <- list(c(1:1050),
                       c(1:322),
                       c(1:2))

yrep2 <- as.data.frame.table(yrep) %>%
  rename("sample" = "Var1",
         'Nest_ID' = "Var2",
         "interval" = "Var3",
         'yrep' = "Freq") %>%
  unite(col = "ID_interval",
        c("Nest_ID", "interval"),
        sep = "_")

conf_df <- yrep2 %>%
  left_join(resp, by = "ID_interval")

conf_summary <- conf_df %>%
  mutate(category = case_when((yrep == 1 & resp == 1) ~ "True Positive",
                              (yrep == 0 & resp == 0) ~ "True Negative",
                              (yrep == 0 & resp == 1) ~ "False Negative",
                              (yrep == 1 & resp == 0) ~ "False Positive",
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
                           'confusion_mod2daily_subset.RDS'))

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

# Summaries ---------------------------------------------------------------

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
             'Model2daily_subset_posterior_samples.RDS'))

summary <- summary(mod$samples)

saveRDS(summary, 
        here('data_outputs',
             '01_whwonests',
             '04_posterior_summaries',
             'Model2daily_subset_posterior_summary.RDS'))

