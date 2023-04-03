# Graphical posterior predictive check
# April 1, 2022
# Ana Miller-ter Kuile

# this is a hack of the bayesplot functionality to generate posterior
# predictive check graphs - specifically to assess - is the model family and link
# function I've selected appropriate for the data I have, or do I need to consider
# a different link or distribution (e.g. logit vs. cloglog link for binomial data; 
# poisson vs. negative binomial distribution for count data)

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "coda", "bayesplot",
                  "jagsUI",
                  "reshape2", "BayesPostEst")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

source(here("code",
            "functions",
            "plot_functions.R"))



# Load GOF model runs -----------------------------------------------------

logit <- readRDS(here('monsoon',
                      'outputs',
                      'model_logit_GOF_2_27_23.RDS'))

cloglog <- readRDS(here('monsoon',
                        'outputs',
                        'model_cloglog_GOF_2_27_23.RDS'))


# Load data ---------------------------------------------------------------

#and we also need our original y data
data <- read.csv(here("data",
                      "Nest_survival_data.csv"))

# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- data %>%
  #set fate to 1-=
  mutate(Fate_class = case_when(Fate_cat == "success" ~ 1,
                                Fate_cat == "failure" ~ 0,
                                TRUE ~ NA_real_)) %>%
  ungroup() %>%
  dplyr::select(Fate_class, Nest_ID) %>%
  #make this type "observed"
  mutate(type = "Observed") %>%
  group_by(Nest_ID) %>%
  mutate(Interval_num = 1:n()) %>%
  ungroup() %>%
  mutate(Nest_ID =  as.numeric(as.factor(Nest_ID)))


# LOGIT -------------------------------------------------------------------

# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps_logit <- logit$sims.list$yrep

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep_logit <- melt(yreps_logit) %>%
  rename("Iteration" = "Var1",
         "Nest_ID" = "Var2",
         "Interval_num" = "Var3",
         "Fate_class" = "value") %>%
  mutate(type = "Simulated")


# Graph observed versus simulated -----------------------------------------


#posterior predictive check graphical observation
ggplot() +
  #graph the simulated data
  geom_density(data = yrep_logit, aes(x = Fate_class, group = Iteration, fill = type), 
               alpha = 0.2) +
  geom_density(data = y, aes(x = Fate_class, fill = type), alpha = 0.5)

#more ones than in our dataset are predicted by this model... hmmm


# CLOGLOG -----------------------------------------------------------------


# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps_cloglog <- cloglog$sims.list$yrep

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep_cloglog <- melt(yreps_cloglog) %>%
  rename("Iteration" = "Var1",
         "Nest_ID" = "Var2",
         "Interval_num" = "Var3",
         "Fate_class" = "value") %>%
  mutate(type = "Simulated")


# Graph observed versus simulated -----------------------------------------


#posterior predictive check graphical observation
ggplot() +
  #graph the simulated data
  geom_density(data = yrep_cloglog, aes(x = Fate_class, group = Iteration, fill = type), 
               alpha = 0.2) +
  geom_density(data = y, aes(x = Fate_class, fill = type), alpha = 0.5)

#more ones than in our dataset are predicted by this model... hmmm

# How often does observed match predicted? --------------------------------

y3 <- y %>%
  rename("observed" = "Fate_class") %>%
  dplyr::select(-type)

# LOGIT -------------------------------------------------------------------

yrepl2 <- yrep_logit %>%
  rename("predicted" = "Fate_class") %>%
  dplyr::select(-Iteration, -type) 

ys_logit <- y3 %>%
  left_join(yrepl2, by = c("Nest_ID", "Interval_num"))

match_0l <- ys_logit %>%
  filter((observed == 0 & observed == predicted)) %>%
  summarise(match_0 = n()) %>%
  as_vector()

mismatch_0l <- ys_logit %>%
  filter((observed == 0 & observed!= predicted)) %>%
  summarise(mismatch_0 = n()) %>%
  as_vector()

match_1l <- ys_logit %>%
  filter((observed == 1 & observed == predicted)) %>%
  summarise(match_1 = n()) %>%
  as_vector()

mismatch_1l <- ys_logit %>%
  filter((observed == 1 & observed != predicted)) %>%
  summarise(mismatch_1 = n()) %>%
  as_vector()

df_logit <- as.data.frame(cbind(type = type <- c("Alive", "Dead"),
match = match <- c(match_1l, match_0l),
mismatch = mismatch <- c(mismatch_1l, mismatch_0l))) %>%
  remove_rownames() %>%
  column_to_rownames(var = "type") %>%
  mutate_if(is.character, as.numeric)

test <- chisq.test(df_logit)

df_logit %>%
  rowwise() %>%
  rownames_to_column(var = "status") %>%
  mutate(accuracy = match/(match+mismatch))

t_logit <- ys_logit %>%
  group_by(Nest_ID, Interval_num, observed) %>%
  summarise(match = sum(observed == predicted, na.rm = T),
            mismatch = sum(observed != predicted, na.rm = T),
            total = match + mismatch,
            perc_match = match/total,
            perc_mismatch = mismatch/total)

t_logit %>%
  mutate(observed_cat = case_when(observed == 1 ~ "alive",
                                  observed == 0 ~ "dead",
                                  TRUE ~ NA_character_)) %>%
  group_by(observed_cat) %>%
  summarise(mean = mean(perc_match, na.rm = T),
            sd = sd(perc_match, na.rm = T),
            total = n(),
            se = sd/sqrt(total)) %>%
  ggplot(aes(x = observed_cat, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2)


# CLOGLOG -----------------------------------------------------------------

yrepc2 <- yrep_cloglog %>%
  rename("predicted" = "Fate_class") %>%
  dplyr::select(-Iteration, -type) 

ys_cloglog <- y3 %>%
  left_join(yrepc2, by = c("Nest_ID", "Interval_num"))

match_0c <- ys_cloglog %>%
  filter((observed == 0 & observed == predicted)) %>%
  summarise(match_0 = n()) %>%
  as_vector()

mismatch_0c <- ys_cloglog %>%
  filter((observed == 0 & observed!= predicted)) %>%
  summarise(mismatch_0 = n()) %>%
  as_vector()

match_1c <- ys_cloglog %>%
  filter((observed == 1 & observed == predicted)) %>%
  summarise(match_1 = n()) %>%
  as_vector()

mismatch_1c <- ys_cloglog %>%
  filter((observed == 1 & observed != predicted)) %>%
  summarise(mismatch_1 = n()) %>%
  as_vector()

df_cloglog <- as.data.frame(cbind(type = type <- c("Alive", "Dead"),
                                match = match <- c(match_1c, match_0c),
                                mismatch = mismatch <- c(mismatch_1c, mismatch_0c))) %>%
  remove_rownames() %>%
  column_to_rownames(var = "type") %>%
  mutate_if(is.character, as.numeric)

test <- chisq.test(df_cloglog)

df_cloglog %>%
  rowwise() %>%
  rownames_to_column(var = "status") %>%
  mutate(accuracy = match/(match+mismatch))

t_cloglog <- ys_cloglog %>%
  group_by(Nest_ID, Interval_num, observed) %>%
  summarise(match = sum(observed == predicted, na.rm = T),
            mismatch = sum(observed != predicted, na.rm = T),
            total = match + mismatch,
            perc_match = match/total,
            perc_mismatch = mismatch/total)

t_cloglog %>%
  mutate(observed_cat = case_when(observed == 1 ~ "alive",
                                  observed == 0 ~ "dead",
                                  TRUE ~ NA_character_)) %>%
  group_by(observed_cat) %>%
  summarise(mean = mean(perc_match, na.rm = T),
            sd = sd(perc_match, na.rm = T),
            total = n(),
            se = sd/sqrt(total)) %>%
  ggplot(aes(x = observed_cat, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2)

# Residuals explorations --------------------------------------------------

#load the formatted data for the JAGS model
source(here("code", 
            "02_nest_survival",
            "01_data_prep.R"))

Nest_ID <- rownames(y)

resid <- model$mean$resid %>%
  as.data.frame() %>%
  cbind(Nest_ID) %>%
  pivot_longer(cols = V1:V15,
               names_to = "interval",
               values_to = "residual") %>%
  mutate(interval = str_sub(interval, 2, length(interval))) %>%
  filter(!is.na(residual)) %>%
  mutate(interval = as.integer(interval)) %>%
  left_join(nests1, by = c("Nest_ID", "interval"))

m1 <- lm(residual ~ Project_ID,
         data = resid)
summary(m1)

siter <- paste("R^2 == 0.001")
(resid_site <- ggplot(resid, aes(x = Project_ID, y = residual)) +
    geom_boxplot() +
    #geom_jitter(height = 0) +
    geom_hline(yintercept = 0, linetype = 2) +
    annotate(geom = "text", 
             x = 3, 
             y = 1,
             label = siter,
             parse = T))

#export this
ggsave(here("figures", 
            "supplement",
            "site_residuals.png"),
       plot = resid_site,
       units = "in",
       height = 6,
       width = 8)

# Patterns in transect-level b0s ------------------------------------------

names <- nests1 %>%
  distinct(Transect_ID, Project_ID) %>%
  mutate(transect = 1:n())

b0 <- param_mod$sims.list$b0
b0.transect <- param_mod$sims.list$b0.transect

b0.transect.t <- b0.transect - b0

b0.transect.df <- as.data.frame(b0.transect.t) %>%
  pivot_longer(cols = V1:V59,
               names_to  = "Transect.num",
               values_to = "b0") %>%
  mutate(Transect.num = 
           str_sub(Transect.num, 2, length(Transect.num))) %>%
  mutate(Transect.num = as.integer(Transect.num)) %>%
  left_join(names, by = c("Transect.num" = "transect")) %>%
  group_by(Transect.num, Project_ID, Transect_ID) %>%
  summarise(median = median(b0, na.rm = T),
            median.LCI = quantile(b0, prob = 0.025, type = 8, na.rm = T),
            median.UCI = quantile(b0, prob = 0.975, type = 8, na.rm = T))  

ggplot(b0.transect.df, aes(x = Transect_ID, y = median, color = Project_ID)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = median.LCI, ymax = median.UCI), width = 0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# AUC ---------------------------------------------------------------------

#this section of script calculates AUC 

#getting a dataset with names that match model (I don't think
# this is necessary, but did it anyway)
nests2 <- nests1 %>%
  mutate(y = case_when(Fate_cat == "success" ~ 1,
                       Fate_cat == "failure" ~ 0,
                       TRUE ~ NA_real_)) %>%
  mutate(Tmax2 = Tmax^2) %>%
  mutate(Trt_cat = factor(Trt_cat, levels = c("U", "B", "H", "HB"))) %>%
  mutate(Tree_sp = factor(Tree_sp, levels = c("PIPO", "Abies", "POTR5",
                                              "JUOC", "PSME"))) %>%
  mutate(Time_groups = factor(Time_groups, 
                              levels = c("oot", "0-3", "4-9", "10+"))) %>%
  mutate(prevStage = factor(prevStage,
                            levels = c("Ne", "Eg"))) %>%
  rowwise() %>%
  mutate(Age = Julian_end - Init_day) %>%
  mutate(PPT2 = PPT^2,
         Age2 = Age^2) %>%
  rename("StageID" = "prevStage",
         "TreatmentID" = "Trt_cat",
         "NestHt" = "Nest_Ht",
         "SpeciesID" = "Tree_sp",
         "InitDay" = "Init_day",
         "Trees50" = "Trees_50",
         "Trees2550" = "Trees_2550",
         "PercPonderosa" = "pPIPO",
         "ForestCV" = "a1000_areacv2",
         "Contag" = "a1000_contag",
         "OpenNm" = "a1000_np1",
         "LandBu" = "a1000_Bu",
         "LandHa" = "a1000_Ha")

colnames(nests2)
# a matrix of the data with all the variables that 
# go into the model
mod_matrix <- model.matrix(y ~   TreatmentID + 
                             n_tx + Age + Age2 + Time_groups + 
                             NestHt +
                             cosOrientation + SpeciesID +
                             InitDay + Trees50 + Trees2550 + 
                             PercPonderosa +
                             Trees50 * PercPonderosa +
                             Trees2550 * PercPonderosa + 
                             Trees50*Tmax +
                             Trees2550*Tmax + 
                             Tmax +  PPT  + PPT2 + ForestCV + 
                             Contag + OpenNm + LandHa + LandBu +
                             LandHa*LandBu + Tmax*PPT,
                           data = nests2)

#a matrix of the data with all the *data* that will go into
# the model
x_data <- as.matrix(model.frame(y ~    TreatmentID + 
                                  n_tx + Age + Age2 + Time_groups + 
                                  NestHt +
                                  cosOrientation + SpeciesID +
                                  InitDay + Trees50 + Trees2550 + 
                                  PercPonderosa +
                                  Trees50 * PercPonderosa +
                                  Trees2550 * PercPonderosa + 
                                  Trees50*Tmax +
                                  Trees2550*Tmax + 
                                  Tmax +  PPT  + PPT2 + ForestCV + 
                                  Contag + OpenNm + LandHa + LandBu +
                                  LandHa*LandBu + Tmax*PPT,
                                data = nests2))

#converting the model with beta estiates to a mcmc object
as.mcmc <- coda::as.mcmc(param_mod)

#defining the variables that we want to pull out of that for
# the mcmc matrix that matches the model and data matrices
vars <- c('b0',
          'b2TreatmentID', "b3TrtTime",
          "b4SpeciesID", 'b'
)
#getting just the simulated list of variables from the mmcmc:
mcmcout <- as.mcmc$sims.list
#selecting those in ou list of variables for the model
mcmcout <- mcmcout[names(mcmcout) %in% vars] 
#making this into a matrix, making sure that the 
# order of the columns matches that of the data and
# model matrices
mcmc_matrix <- as.data.frame(do.call(cbind, mcmcout)) %>%
  dplyr::select(b0, 
                V3:V5, 
                V7:V9,
                V11:V14,
                V17:V39) %>%
  as.matrix()

# using mcmcRocPrcGen, determining the ROC curve for the 
#model
fit_sum1 <- mcmcRocPrcGen(modelmatrix = mod_matrix,
                          modelframe = x_data,
                          mcmcout = mcmc_matrix,
                          curves = TRUE,
                          fullsims = FALSE)

#getting the AUC stat from that output
fit_sum1$area_under_roc

#plotting that ROC curve
ggplot(data = as.data.frame(fit_sum1$roc_dat), 
       aes(x = x, y = y)) +
  geom_line() + 
  geom_abline(intercept = 0, slope = 1, color = "gray") + 
  labs(title = "ROC curve") + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") +
  annotate(geom = "text", x = 0.75, y = 0.25, label = "AUC = 0.49")


