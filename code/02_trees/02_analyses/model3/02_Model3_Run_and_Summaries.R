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
                  "mcmcplots") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------
# 
# #load the formatted data for the JAGS model
data <- readRDS(here("data_outputs",
                     "02_trees",
                     "03_JAGS_input_data",
                     "mod3_JAGS_input_data.RDS"))

# Parameters to save ------------------------------------------------------


params <- c('b0',
            #covariates
            'b',
            'b1'
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "02_trees",
              "02_analyses",
              "model3",
              "jags",
              "model3.R")

mod <- jagsUI::jags(data = data,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    n.thin = 3,
                    n.burnin = 1000,
                    DIC = TRUE)

mcmcplot(mod$samples)

gelman.diag(mod$samples,
            multivariate = F)

# Update with initials ----------------------------------------------------
# 
# #initials for next run
# #get the MCMC chains
# samples <- mod$samples
# 
# #function to make each chain a dataframe
# df_fun <- function(chain){
#   df <- as.data.frame(chain) %>%
#     rownames_to_column(var = "sample")
#   return(df)
# }
# 
# #use that function on all list elements
# samp_dfs <- lapply(samples, df_fun)
# 
# #make into one dataframe
# samp_df <- bind_rows(samp_dfs, .id = "chain")
# 
# #get values for all parameters from the last iteration of the
# #chain with the lowest deviance
# samp_df2 <- samp_df %>%
#   group_by(chain) %>%
#   #get mean deviance by chain
#   mutate(mean_dev = mean(deviance, na.rm = T)) %>%
#   ungroup() %>%
#   #get only the chain with the minimum average deviance
#   filter(mean_dev == min(mean_dev)) %>%
#   #pull out the final sample from that chain
#   filter(sample == max(sample)) %>%
#   dplyr::select(-chain, -sample,
#                 -deviance, -mean_dev)
# 
# b0 <- samp_df2 %>%
#   dplyr::select(b0) %>%
#   as_vector()
# 
# b1 <- samp_df2 %>%
#   dplyr::select(`b1[2]`) %>%
#   as_vector()
# 
# b1 <- c(NA, b1)
# 
# b <- samp_df2 %>%
#   dplyr::select(`b[2]`:`b[7]`) %>%
#   as_vector()
# 
# b <- c(NA, b)
# # 
# inits <- list(list(b0 = b0,
#                    b1= b1,
#                    b = b),
#               list(b0 = b0 +0.1,
#                    b = b+0.1,
#                    b1 = b1+0.1),
#               list(b0 = b0-0.1,
#                    b = b-0.1,
#                    b1 = b1-0.1))
# 
# mod2 <- jagsUI::jags(data = data,
#                     inits = inits,
#                     model.file = model,
#                     parameters.to.save = params,
#                     parallel = TRUE,
#                     n.chains = 3,
#                     n.iter = 10000,
#                     n.thin = 4,
#                     n.burnin = 2000,
#                     DIC = TRUE)
# 
# mcmcplot(mod2$samples)
# 
# gelman.diag(mod2$samples,
#             multivariate = F)
# 
# 
# # Update again with inits -------------------------------------------------
# 
# #initials for next run
# #get the MCMC chains
# samples2 <- mod2$samples
# 
# #use that function on all list elements
# samp_dfs <- lapply(samples2, df_fun)
# 
# #make into one dataframe
# samp_df <- bind_rows(samp_dfs, .id = "chain")
# 
# #get values for all parameters from the last iteration of the
# #chain with the lowest deviance
# samp_df2 <- samp_df %>%
#   group_by(chain) %>%
#   #get mean deviance by chain
#   mutate(mean_dev = mean(deviance, na.rm = T)) %>%
#   ungroup() %>%
#   #get only the chain with the minimum average deviance
#   filter(mean_dev == min(mean_dev)) %>%
#   #pull out the final sample from that chain
#   filter(sample == max(sample)) %>%
#   dplyr::select(-chain, -sample,
#                 -deviance, -mean_dev)
# 
# b0 <- samp_df2 %>%
#   dplyr::select(b0) %>%
#   as_vector()
# 
# b1 <- samp_df2 %>%
#   dplyr::select(`b1[2]`) %>%
#   as_vector()
# 
# b1 <- c(NA, b1)
# 
# b <- samp_df2 %>%
#   dplyr::select(`b[2]`:`b[7]`) %>%
#   as_vector()
# 
# b <- c(NA, b)
# # 
# inits <- list(list(b0 = b0,
#                    b1= b1,
#                    b = b),
#               list(b0 = b0 +0.1,
#                    b = b+0.1,
#                    b1 = b1+0.1),
#               list(b0 = b0-0.1,
#                    b = b-0.1,
#                    b1 = b1-0.1))
# 
# mod3 <- jagsUI::jags(data = data,
#                      inits = inits,
#                      model.file = model,
#                      parameters.to.save = params,
#                      parallel = TRUE,
#                      n.chains = 3,
#                      n.iter = 10000,
#                      n.thin = 4,
#                      n.burnin = 2000,
#                      DIC = TRUE)
# 
# mcmcplot(mod3$samples)
# 
# gelman.diag(mod3$samples,
#             multivariate = F)

Rhat <- mod$Rhat

saveRDS(Rhat,
        here('data_outputs',
        '02_trees',
        '04_posterior_summaries',
        'model3',
        'model3_Rhat.RDS'))


# Update for goodness of fit ----------------------------------------------

parms <- c("yrep_1", "yrep_2", 'resid_1', "resid_2", 'p1', "p2")

mod.update <- update(mod,
                     parameters.to.save = parms,
                     n.iter = 350,
                     codaOnly = TRUE)

saveRDS(mod.update,
        here('data_outputs',
             '02_trees',
             '04_posterior_summaries',
             'model3',
             'GOF_model3.RDS'))
        

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
             '02_trees',
             '04_posterior_summaries',
             'model3',
             'samples_model3.RDS'))

summary <- summary(mod$samples)

saveRDS(summary,
        here('data_outputs',
             '02_trees',
             '04_posterior_summaries',
             'model3',
             'summary_model3.RDS'))

