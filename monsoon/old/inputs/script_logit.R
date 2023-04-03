# Running the nest survival model
# Ana Miller-ter Kuile
# November 4, 2021

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
data <- readRDS("/scratch/atm234/nest_survival/inputs/survival_JAGS_input_data.RDS")

# Compile data ------------------------------------------------------------


data_list <- list(#overall values for likelihood loops
                  n.nests = data$n.nests,
                  n.t = data$n.t,
                  n.forests = data$n.forests,
                  #Random variables
                  Nest.num = data$Nest.num,
                  Transect.num = data$Transect.num,
                  Year.num = data$Year.num,
                  Forest.num = data$Forest.num,
                  #Interval covariates
                  StageID = data$Stage,
                  #Treatment covariate
                  TreatmentID = data$TreatmentID,
                  #Nest-level covariates
                  NestHt = data$NestHt,
                  cosOrientation = data$cosOrientation,
                  InitDay = data$InitDay,
                  SpeciesID = data$SpeciesID,
                  #local-level
                  Trees50 = data$Trees50,
                  Trees2550 = data$Trees2550,
                  PercPonderosa = data$PercPonderosa,
                  #landscape-level
                  Tmax = data$Tmax,
                  LandHa = data$LandHa,
                  LandBu = data$LandBu,
                  ForestCV = data$ForestCV,
                  Contag = data$Contag,
                  OpenNm = data$OpenNm,
                  #Data
                  y = data$y,
                  t = data$t,
                  #numbers for prior distribution loops
                  n.transects = data$n.transects,
                  n.years = data$n.years,
                  n.trt = data$n.trt,
                  n.stages = data$n.stages,
                  n.species = data$n.species)


# Parameters to save ------------------------------------------------------



params <- c(#Environmental covariates
            'b1StageID',
            'b2TreatmentID',
            'b3SpeciesID',
            'b',
            #Random covariate betas
            'b0',
            'b0.transect',
            'b0.nest',
            'b0.year',
            'sig.transect',
            'sig.year',
            'sig.nest',
            #p-value objects
            'z.b1',
            'z.b2',
            'z.b3',
            'z')
          
          

# INits -------------------------------------------------------------------


# JAGS model --------------------------------------------------------------

egg_mod <- jagsUI::jags(data = data_list,
                        #inits = inits,
                        inits = NULL,
                        model.file = "/scratch/atm234/nest_survival/inputs/model_logit.R",
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.burnin = 1000,
                        n.iter = 294000,
                        #n.thin = 10,
                        DIC = TRUE)

#save as an R data object
saveRDS(egg_mod, 
        file = "/scratch/atm234/nest_survival/outputs/model_logit_MCMC_2_8_23.RDS")

Sys.time()
