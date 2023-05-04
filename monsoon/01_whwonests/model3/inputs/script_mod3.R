#Monsoon script - Model3: Conditional probability interval model
# Ana Miller-ter Kuile
# April 5, 2023

#this script runs the custom interval probability model

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
data <- readRDS("/scratch/atm234/survival_models/nests/model3/inputs/mod3_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(#Data count variables
                  n.nests1 = data$n.nests1,
                  n.nests = data$n.nests,
                  n.t = data$n.t, 
                  n.years = data$n.years,
                  n.transects = data$n.transects,
                  n.trt = data$n.trt, 
                  n.species = data$n.species,
                  n.stages = data$n.stages,
                  n.forests = data$n.forests,
                  #Random effects IDs
                  Year.num = data$Year.num,
                  Transect.num = data$Transect.num,
                  Forest.num = data$Forest.num,
                  #missing data
                  Forest.ID = data$Forest.ID,
                  #Treatment covariate
                  TreatmentID = data$TreatmentID, 
                  #nest tree species
                  SpeciesID = data$SpeciesID, 
                  #Interval stage covariate
                  StageID = data$Stage,
                  #Nest-level covariates
                  NestHt = data$NestHt, 
                  cosOrientation = data$cosOrientation,
                  InitDay = data$InitDay, 
                  #Local-level covariates
                  Trees50 = data$Trees50,
                  Trees2550 = data$Trees2550, 
                  PercPonderosa = data$PercPonderosa,
                  #climate covariates
                  Tmax = data$Tmax,
                  PPT = data$PPT,
                  #landscape-scale covariates
                  ForestCV = data$ForestCV,
                  Contag = data$Contag,
                  OpenNm = data$OpenNm,
                  LandHa = data$LandHa,
                  LandBu = data$LandBu,
                  #dataset
                  y = data$y,
                  #duplicated dataset for ifelse statement
                  y2 = data$y,
                  #interval lengths
                  t = data$t)
# Parameters to save ------------------------------------------------------
params <- c(
            #Random covariate betas
            'b0.forest',
            'b0.transect',
            'b0.year',
            'b0',
            #Variance/precision
            'sig.forest',
            'sig.transect',
            'sig.year',
            #covariates
            'b1TreatmentID',
            'b2SpeciesID',
            "b3StageID",
            'b',
            'z.b1',
            'z.b2',
            'z.b3',
            'z'
          )

      

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = NULL,
                        model.file = "/scratch/atm234/survival_models/nests/model3/inputs/model3.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.burnin = 1000,
                        n.iter = 100000,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file = "/scratch/atm234/survival_models/nests/model3/outputs/model3_JAGS_model.RDS")

Sys.time()


