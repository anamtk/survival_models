# Run high var JAGS models
# Ana Miller-ter Kuile
# May 17, 2023

#script to run models for simulated data with high predictor variation

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  "jagsUI",
                  "coda",
                  "mcmcplots")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data1 <- readRDS(here("data_outputs",
                    "03_simulated",
                    "03_jags_input_data",
                    "Model1high_JAGS_data.RDS"))


data2 <- readRDS(here("data_outputs",
                    "03_simulated",
                    "03_jags_input_data",
                    "Model2high_JAGS_data.RDS"))

data3 <- readRDS(here("data_outputs",
                    "03_simulated",
                    "03_jags_input_data",
                    "Model3high_JAGS_data.RDS"))

# Set up for jags models --------------------------------------------------

parms <- c("b1", "b0")

# Run model 1 -------------------------------------------------------------

model1 <- here("code",
               "03_simulated_datasets",
               "02_analyses",
               "00_jags",
               "model1.R")

Sys.time()
mod1 <- jagsUI::jags(data = data1,
                     inits = NULL,
                     model.file = model1,
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter =15000,
                     n.thin = 3,
                     DIC = TRUE)
Sys.time()

gel1 <- gelman.diag(mod1$samples)
as.data.frame(gel1$psrf) %>%
  ggplot(aes(x = `Point est.`)) +
  geom_histogram() +
  geom_vline(xintercept = 1.1) +
  geom_vline(xintercept = 1.2)

#mcmcplot(mod1$samples)
# 
sum1 <- summary(mod1$samples)
# 
saveRDS(sum1, here('data_outputs',
                   '03_simulated',
                   '04_posterior_summaries',
                   'mod1high_summary.RDS'))
 
parms1 <- c("yrep", 'p')

# update number of iterations to always
#get 1050 posterior samples
mod1update <- update(mod1,
                     parameters.to.save = parms1,
                     n.iter = 1051)
# 
saveRDS(mod1update, here('data_outputs',
                         '03_simulated',
                         '04_posterior_summaries',
                         'mod1high_GOFmodel.RDS'))

# Run model 2 -------------------------------------------------------------

model2 <- here("code",
               "03_simulated_datasets",
               "02_analyses",
               "00_jags",
               "model2.R")

Sys.time() 
mod2 <- jagsUI::jags(data = data2,
                     inits = NULL,
                     model.file = model2,
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 25000,
                     n.thin = 10,
                     DIC = TRUE)


gelman.diag(mod2$samples)
gel2s <- gelman.diag(mod2$samples)
as.data.frame(gel2s$psrf) %>%
  ggplot(aes(x = `Point est.`)) +
  geom_histogram() +
  geom_vline(xintercept = 1.1) +
  geom_vline(xintercept = 1.2)

sum2 <- summary(mod2$samples)

saveRDS(sum2, here('data_outputs',
                   '03_simulated',
                   '04_posterior_summaries',
                   'mod2high_summary.RDS'))

parms2 <- c("y.repkeep", 'p.intkeep', 'p.int', 'yrep')

mod2update <- update(mod2,
                     parameters.to.save = parms2,
                     n.iter = 3500)
# 
saveRDS(mod2update, here('data_outputs',
                         '03_simulated',
                         '04_posterior_summaries',
                        'mod2high_GOFmodel.RDS'))

# Run model 3 -------------------------------------------------------------

model3 <- here("code",
               "03_simulated_datasets",
               "02_analyses",
               "00_jags",
               "model3_April2025_daily2.R")

mod3 <- jagsUI::jags(data = data3,
                     inits = NULL,
                     model.file = model3,
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.burnin = 2000,
                     n.iter = 100000,
                     n.thin = 10,
                     DIC = TRUE)

gelman.diag(mod3$samples)
gel3s <- gelman.diag(mod3$samples)
as.data.frame(gel3s$psrf) %>%
  ggplot(aes(x = `Point est.`)) +
  geom_histogram() +
  geom_vline(xintercept = 1.1)+
  geom_vline(xintercept = 1.2)

#initials for next run
#get the MCMC chains
# samples <- mod3$samples
# 
# #function to make each chain a dataframe
# df_fun <- function(chain){
#   df <- as.data.frame(chain) %>%
#     rownames_to_column(var = "iteration")
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
#   #pull out the final iteration from that chain
#   filter(iteration == max(iteration)) %>%
#   dplyr::select(-chain, -iteration,
#                 -deviance, -mean_dev) 
# 
# b1 <- samp_df2 %>%
#   dplyr::select(`b1[1]`:`b1[100]`) %>%
#   as_vector()
# 
# b0 <- samp_df2 %>%
#   dplyr::select(`b0[1]`:`b0[100]`) %>%
#   as_vector()
# 
# inits3 <- list(list(b0 = b0,
#                         b1= b1),
#                    list(b0 = b0 +0.05,
#                         b1 = b1+0.05),
#                    list(b0 = b0+0.08,
#                         b1 = b1+0.08))
# 
# mod3_2 <- jagsUI::jags(data = data3,
#                        inits = inits3,
#                        model.file = model3,
#                        parameters.to.save = parms,
#                        parallel = TRUE,
#                        n.chains = 3,
#                        n.burnin = 2000,
#                        n.iter = 100000,
#                        n.thin = 10,
#                        DIC = TRUE)
# 
# gelman.diag(mod3_2$samples)
# gel32s <- gelman.diag(mod3_2$samples)
# as.data.frame(gel32s$psrf) %>%
#   ggplot(aes(x = `Point est.`)) +
#   geom_histogram() +
#   geom_vline(xintercept = 1.1)+
#   geom_vline(xintercept = 1.2)

sum3 <- summary(mod3$summary)

saveRDS(sum3, here('data_outputs',
                   '03_simulated',
                   '04_posterior_summaries',
                   'mod3high_summary.RDS'))
# 
parms3 <- c("yrep_1", "yrep_2",'pi1', "pi2")

mod3update <- update(mod3,
                     parameters.to.save = parms3,
                     n.iter = 3500)
# 
saveRDS(mod3update, here('data_outputs',
                         '03_simulated',
                         '04_posterior_summaries',
                         'mod3high_GOFmodel.RDS'))

