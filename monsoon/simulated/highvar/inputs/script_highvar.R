# Run high var JAGS models
# Ana Miller-ter Kuile
# May 17, 2023

#script to run models for simulated data with high predictor variation

# Load packages -----------------------------------------------------------

package.list <- c("jagsUI",
                  "coda",
                  "mcmcplots")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data1 <- readRDS("/scratch/atm234/survival_models/sim/highvar/inputs/Model1high_JAGS_data.RDS")


data2 <- readRDS("/scratch/atm234/survival_models/sim/highvar/inputs/Model2high_JAGS_data.RDS")


data3 <- readRDS("/scratch/atm234/survival_models/sim/highvar/inputs/Model3high_JAGS_data.RDS")


# Set up for jags models --------------------------------------------------

parms <- c("b1", "b0")

# Run model 1 -------------------------------------------------------------

Sys.time()
print("startmod1")

data1_list <- list(n.datasets = data1$n.datasets,
                   n.indiv = data1$n.indiv,
                   y = data1$y,
                   x = data1$x,
                   t = data1$t)

mod1 <- jagsUI::jags(data = data1_list,
                     inits = NULL,
                     model.file = "/scratch/atm234/survival_models/sim/models/model1.R",
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 4000,
                     DIC = TRUE)

mcmcplot(mod1$samples,
         dir = '/scratch/atm234/survival_models/sim/highvar/outputs/mcmcplots/model1')

sum1 <- summary(mod1$samples)

saveRDS(sum1, 
        '/scratch/atm234/survival_models/sim/highvar/outputs/mod1_highvar_sum.RDS')

parms1 <- c( "p")

mod1update <- update(mod1,
                     parameters.to.save = parms1,
                     n.iter = 350)

saveRDS(mod1update, 
        '/scratch/atm234/survival_models/sim/highvar/outputs/mod1_GOF.RDS')

Sys.time()
print("endmod1")
# Run model 2 -------------------------------------------------------------

Sys.time()
print("startmod2")

data2_list <- list(n.datasets = data2$n.datasets,
                   n.indiv = data2$n.indiv,
                   n.t = data2$n.t,
                   y = data2$y,
                   x = data2$x,
                   t = data2$t)


mod2 <- jagsUI::jags(data = data2_list,
                     inits = NULL,
                     model.file = "/scratch/atm234/survival_models/sim/models/model2.R",
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 4000,
                     DIC = TRUE)

mcmcplot(mod2$samples,
         dir = '/scratch/atm234/survival_models/sim/highvar/outputs/mcmcplots/model2')

sum2 <- summary(mod2$samples)

saveRDS(sum2, 
        '/scratch/atm234/survival_models/sim/highvar/outputs/mod2_highvar_sum.RDS')

parms2 <- c("p.intkeep", "p.int")

mod2update <- update(mod2,
                     parameters.to.save = parms2,
                     n.iter = 350)

saveRDS(mod2update,
        '/scratch/atm234/survival_models/sim/highvar/outputs/mod2_GOF.RDS')

Sys.time()
print("endmod2")
# Run model 3 -------------------------------------------------------------

Sys.time()
print("startmod3")

data3_list <- list(n.datasets = data3$n.datasets,
                   n.indiv = data3$n.indiv,
                   n.indiv1 = data3$n.indiv1,
                   n.t = data3$n.t,
                   y = data3$y,
                   x = data3$x,
                   t = data3$t)

mod3 <- jagsUI::jags(data = data3_list,
                     inits = NULL,
                     model.file = "/scratch/atm234/survival_models/sim/models/model3.R",
                     parameters.to.save = parms,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 4000,
                     DIC = TRUE)

mcmcplot(mod3$samples,
         dir = '/scratch/atm234/survival_models/sim/highvar/outputs/mcmcplots/model3')

sum3 <- summary(mod3$samples)

saveRDS(sum3, '/scratch/atm234/survival_models/sim/highvar/outputs/mod3_highvar_sum.RDS')

parms3 <- c("p1", "p2")

mod3update <- update(mod3,
                     parameters.to.save = parms3,
                     n.iter = 350)

saveRDS(mod3update, '/scratch/atm234/survival_models/sim/highvar/outputs/mod3_GOF.RDS')

Sys.time()
print("endmod3")
