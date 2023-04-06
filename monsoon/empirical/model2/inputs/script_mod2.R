#Monsoon script - Model2: normal daily survival model
# Ana Miller-ter Kuile
# April 5, 2023

#this script runs the original interval daily survival model

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages
package.list <- c("jagsUI", "coda") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/survival_models/model2/inputs/mod2_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(#overall values for likelihood loops
                  #Data count variables
                  n.nests = data$n.nests,
                  n.t = data$n.t, 
                  n.years = data$n.years,
                  n.transects = data$n.transects,
                  n.trt = data$n.trt, 
                  n.species = data$n.species,
                  n.stages = data$n.stages,
                  n.forests = data$n.forests,
                  #Random effects IDs
                  Nest.num = data$Nest.num,
                  Year.num = data$Year.num,
                  Transect.num = data$Transect.num,
                  Forest.num = data$Forest.num,
                  #Interval covariate
                  StageID = data$Stage,
                  Age = data$Age,
                  #Treatment covariate
                  TreatmentID = data$TreatmentID, 
                  #Nest-level covariates
                  NestHt = data$NestHt, 
                  cosOrientation = data$cosOrientation,
                  InitDay = data$InitDay, 
                  SpeciesID = data$SpeciesID, 
                  #Local-level covariates
                  Trees50 = data$Trees50,
                  Trees2550 = data$Trees2550, 
                  PercPonderosa = data$PercPonderosa,
                  #landscape-scale covariates
                  Tmax = data$Tmax,
                  PPT = data$PPT,
                  ForestCV = data$ForestCV,
                  Contag = data$Contag,
                  OpenNm = data$OpenNm,
                  LandHa = data$LandHa,
                  LandBu = data$LandBu,
                  #dataset
                  y = data$y,
                  #interval lengths
                  t = data$t)

# Parameters to save ------------------------------------------------------
params <- c(
            #Random covariate betas
            'b0.nest',
            'b0.transect',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.nest',
            'sig.transect',
            'sig.year',
            'b',
            'b1StageID',
            'b2TreatmentID',
            'b3SpeciesID'
          )

# INits -------------------------------------------------------------------

inits <- readRDS("/scratch/atm234/survival_models/model2/inputs/model2_inits.RDS")

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = inits,
                        #inits = NULL,
                        model.file = "/scratch/atm234/survival_models/model2/inputs/model2.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.burnin = 1000,
                        n.iter = 47000,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file = "/scratch/atm234/survival_models/model2/outputs/model2_JAGS_model.RDS")

Sys.time()


