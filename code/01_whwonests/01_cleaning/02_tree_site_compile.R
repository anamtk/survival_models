# Point-level tree and snag data
# Ana Miller-ter Kuile
# October 29, 2021

# this script combines the tree and snag data into
# summary datasets by survey point

# Steps in this script:
# 1. Figure out the total area of each transect and convert counts to 
# densities (MIS transects are wonky)

# 2. Convert Stage5 large tree DBH (actually Diameter at base) to DBH and 
 # remove those that would not be part of survey
## from Todd's code # PIPO (DBH): y = 0.880x - 3.807, R2 = 0.95, n = 328 (use also for PICO)
## # PSME/PSMEG (DBH): y = 0.853x - 0.996, R2 = 0.94, n = 294 (also use for ABCO and ABGR)

# 3. Remove stage 5 trees from "post" and "never_treated" sites

# 4. Subset snags > 23 (or 25?) cm

# 5. Create a final dataset of survey point and year, total small, med, 
# large and snag densities

# 6. Determine %PIPO for each plot (all trees > 25cm)

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "readxl")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#custom functions 
source(here("code", 
            "functions", 
            "tidy_functions.R"))

# Load datasets -----------------------------------------------------------

trees <- read_xlsx(here("data_raw",
                        "01_whwonests",
                        "veg_data",
                        "Veg02_Live_trees_over_25.xlsx"))

#Info on plot sizes
meta <- read_xlsx(here("data_raw",
                       "01_whwonests",
                       "veg_data",
                       "metadata",
                       "Veg00_Habitat_header_info.xlsx"))

# Determine plot sizes ----------------------------------------------------

#from the metadata for tree plots, subset only EM sites
meta2 <- meta %>%
  filter(str_detect(Point_ID, "EMMAOR|EMFWOR|EMPAID")) #select EM sites only

#want three columns:
# sm25plot default: 2m x 200 m (0.04ha), special cases 2m x 100 (0.02ha)
# md25plot default: 6m x 200m (0.12 ha), special cases 6m x 100 (0.06 ha)
# lg50plot (includes snags) default: 20m x 200 (0.4 ha), special cases 20 x 100 (0.2 ha)

# From Todd: EMPAID transects with different sized plots
# ALL OF EMPAID_SF 
# EMPAID_1WS points 02-06
# ALL OF EMPAID_BR 
# ALL OF EMPAID_CU
# Get a vector of the points that are in the smaller size plot category
points <- meta2 %>%
  #select all the points that are a different size
  filter(str_detect(Point_ID, "EMPAID_SF|EMPAID_BR|EMPAID_CU|EMPAID_1WS")) %>%
  #get distinct points
  distinct(Point_ID) %>%
  #remove points that are not the smaller size 
  filter(!Point_ID %in% c("EMPAID_1WS07", "EMPAID_1WS08", 
                          "EMPAID_1WS09", "EMPAID_1WS10")) %>%
  #make this a vector for the next step of conditional formatting
  as_vector()

#Create plot sizes with a conditional set of criteria
meta3 <- meta2 %>%
  #small trees - when in the points above, make them smaller, otherwise bigger
  mutate(sm25plot = case_when(Point_ID %in% points ~ 0.02, 
                              TRUE ~ 0.04),
         #medium trees - do the same thing for small/large distinction
         md25plot = case_when(Point_ID %in% points ~ 0.06,
                              TRUE ~ 0.12),
         # large trees and all snags - small/large distinctoin
         lg50plot = case_when(Point_ID %in% points ~ 0.2,
                              TRUE ~ 0.4)) %>%
  #from this metadata, select the variables of interest for combining below
  dplyr::select(Measurement_ID, Point_ID, sm25plot, md25plot, lg50plot,
                Trt)


# Stage 5 DBH convserion --------------------------------------------------

# 2. Convert Stage5 large tree DBH (actually Diameter at base) to DBH
## from Todd's code # PIPO (DBH): y = 0.880x - 3.807, R2 = 0.95,
#n = 328 (use also for PICO)
## # PSME/PSMEG (DBH): y = 0.853x - 0.996, R2 = 0.94, 
#n = 294 (also use for ABCO and ABGR)

#eventually, will want to delete any stage 5 that were measured but shouldn't
# be included, which is 
# 3 m for 25 - 49, and 10 m from center for >=50

#with the tree dataset
big_trees <- trees %>%
  #select only areas of effectiveness monitoring
  filter(str_detect(Measurement_ID, "EMPAID|EMFWOR|EMMAOR")) %>%
  #make the distance measured in the comments from center into a number
  mutate(dist = numextract(Comments)) %>% 
  mutate(dist = as.numeric(dist)) %>% #convert to numeric
  #convert stump measurements to DBH given eqn from Todd's code
  mutate(DBH2 = case_when(Struc_class == 5 & 
                           Tree_species %in% c("PIPO", 
                                               "PICO") ~ 0.88*DBH - 3.807, #PInes
                         Struc_class == 5 &
                           Tree_species %in% c("PSME", 
                                               "PSMEG", 
                                               "ABCO", 
                                               "ABGR") ~ 0.853*DBH - 0.996, #others
                         TRUE ~ DBH))# %>% #29360 any not in that structure class just get DBH2 as their DBH

removed_trees <- trees %>%
  #select only areas of effectiveness monitoring
  filter(str_detect(Measurement_ID, "EMPAID|EMFWOR|EMMAOR")) %>%
  #make the distance measured in the comments from center into a number
  mutate(dist = numextract(Comments)) %>% 
  mutate(dist = as.numeric(dist)) %>% #convert to numeric
  #convert stump measurements to DBH given eqn from Todd's code
  mutate(DBH2 = case_when(Struc_class == 5 & 
                            Tree_species %in% c("PIPO", 
                                                "PICO") ~ 0.88*DBH - 3.807, #PInes
                          Struc_class == 5 &
                            Tree_species %in% c("PSME", 
                                                "PSMEG", 
                                                "ABCO", 
                                                "ABGR") ~ 0.853*DBH - 0.996, #others
                          TRUE ~ DBH)) %>% #29360 any not in that structure class just get DBH2 as their DBH
  # remove any now too-small trees in the 3-20m dist
  filter((Struc_class == 5 & dist > 3 & DBH2 < 50)) 

removed_trees2 <- trees %>%
  #select only areas of effectiveness monitoring
  filter(str_detect(Measurement_ID, "EMPAID|EMFWOR|EMMAOR")) %>%
  #make the distance measured in the comments from center into a number
  mutate(dist = numextract(Comments)) %>% 
  mutate(dist = as.numeric(dist)) %>% #convert to numeric
  #convert stump measurements to DBH given eqn from Todd's code
  mutate(DBH2 = case_when(Struc_class == 5 & 
                            Tree_species %in% c("PIPO", 
                                                "PICO") ~ 0.88*DBH - 3.807, #PInes
                          Struc_class == 5 &
                            Tree_species %in% c("PSME", 
                                                "PSMEG", 
                                                "ABCO", 
                                                "ABGR") ~ 0.853*DBH - 0.996, #others
                          TRUE ~ DBH)) %>% #29360 any not in that structure class just get DBH2 as their DBH
  #remove any now too-small trees in the 0-3m dist
  filter((Struc_class == 5 & dist <= 3 & DBH2 < 25 )) %>%#28285
  rbind(removed_trees) 

#anti_join those now too-small trees
big_trees <- big_trees %>%
  anti_join(removed_trees2)

#Note: This point has become somewhat unecessary since Jon
#fixed the database based on the results of this "removed_trees2"
# df, but keepign it in here so we remember we did it

# Get "no tree" points ----------------------------------------------------
#For some reason, the above filter removes the "NONE" category, which
#makes it look like some points weren't surveyed
none <- big_trees %>%
  filter(Tree_species == "NONE") %>%
  distinct(Measurement_ID) %>%
  mutate(Trees2550 = 0,
         Trees50 = 0,
         PIPO2550 = 0,
         PIPO50 = 0) %>%
  dplyr::select(Measurement_ID, Trees2550,
                Trees50, PIPO2550,
                PIPO50) %>%
  ungroup()

# Treatment codes and remove stage 5 in some ------------------------------

# treatment codes
meta3 <- meta3 %>%
  #set treatment period of treated-untreated based on whether "U" or not
  mutate(trt_cat = case_when(Trt != 'U' ~ "treated",
                                TRUE ~ "untreated")) %>%
  #combine all treatment codes here
  mutate(Trt_code = case_when(Trt %in% c("H", "OA", "J", "T",
                                    "HT") ~ "H",
                         Trt %in% c("OHB", "OHTB", "OTHB", "HB", "OBH", "OB",
                                    "B H", "TB", "HTB", "BH") ~ "HB",
                         Trt == "B" ~ "B",
                         Trt == "U" ~ "U",
                         TRUE ~ NA_character_))


#combine with visits dataset
big_trees2 <- big_trees %>%
  left_join(meta3, by = "Measurement_ID") 

#deal with areas that were only visited in post-treatment
# by duplicating them for a pre-treatment period and then deleting
# structure class 5 from the post treatment ones
post_only <- big_trees2 %>%
  group_by(Point_ID) %>%
  #give a category based on the treatment categories indicated
  mutate(group = case_when(all(trt_cat == "treated") ~ "post-only",
                           all(trt_cat == "untreated") ~ "never_treated",
                           TRUE ~ "pre-post")) %>%
  ungroup() %>%
  #select only those surveyed post-treatment
  filter(group == "post-only") %>%
  #give them a new measurement ID with "pre" on it
  mutate(Measurement_ID = paste(Point_ID, "pre", sep = "_")) %>%
  #give them a new treatment period of "untreated"
  mutate(trt_cat = "untreated",
         Trt = "U") %>% # and new treatment code of "U"
  #get rid of the group variable
  dplyr::select(-group)
  
#recombine with the rest of the tree dataset
big_trees3 <- big_trees2 %>%
  #bind the rows
  rbind(post_only) %>% 
  #then delete structure class 5 from anything post treatment or untreated
  filter(!(Struc_class == 5 & Trt != "U"))


## Add new measurement IDs to the metadata frames (visits and meta)
#new poitns to add to meta
points_new <- post_only %>%
  distinct(Measurement_ID, Point_ID) %>%
  mutate(sm25plot = NA, 
         md25plot = NA, 
         lg50plot = NA,
         Trt = "U",
         trt_cat = "untreated",
         Trt_code = "U")

#add new points to meta and fill missing values by point ID 
meta4 <- meta3 %>%
  rbind(points_new) %>%
  group_by(Point_ID) %>%
  fill(sm25plot, md25plot, lg50plot, .direction = "updown") %>%
  ungroup()

# Tree Tallies ------------------------------------------------------------
#how many trees per measurement period?
big_trees3 %>%
  filter(DBH2 < 50) %>%
  filter(Tree_species != c("UNKN")) %>%
  group_by(Measurement_ID) %>%
  tally() %>%
  summarise(mean = mean(n),
            sd = sd(n))

#how many trees per measurement period?
big_trees3 %>%
  filter(DBH2 >= 50) %>%
  filter(Tree_species != c("UNKN")) %>%
  group_by(Measurement_ID) %>%
  tally() %>%
  summarise(mean = mean(n),
            sd = sd(n))

#trees 25-49, measured in the smaller 6 m band
Trees2550 <- big_trees3 %>%
  filter(DBH2 < 50) %>%
  group_by(Measurement_ID) %>%
  tally(name = "Trees2550") %>%
  ungroup()

#only PIPO of that size group
PIPO2550 <- big_trees3 %>%
  filter(DBH < 50) %>%
  filter(Tree_species == "PIPO") %>%
  group_by(Measurement_ID) %>%
  tally(name = "PIPO2550") %>%
  ungroup()

#trees 50 and greater, measured in the larger 20 m band
Trees50 <- big_trees3 %>%
  filter(DBH2 >= 50) %>%
  group_by(Measurement_ID) %>%
  tally(name = "Trees50") %>%
  ungroup()

#PIPO50 and greater
PIPO50 <- big_trees3 %>%
  filter(DBH2 >= 50) %>%
  filter(Tree_species == "PIPO") %>%
  group_by(Measurement_ID) %>%
  tally(name = "PIPO50") %>%
  ungroup()


# Percent PIPO ------------------------------------------------------------

percPIPO <- big_trees3 %>%
  group_by(Measurement_ID) %>%
  mutate(total_trees = n()) %>%
  ungroup() %>%
  group_by(Measurement_ID, total_trees, Tree_species) %>%
  tally() %>%
  filter(Tree_species == "PIPO") %>%
  mutate(pPIPO = n/total_trees) %>%
  ungroup() %>%
  dplyr::select(Measurement_ID, pPIPO)


# Combine all data  -------------------------------------------------------

# add all the DFs together!
tree_tallies <- list(Trees2550, #all the DFs we created of counts
                     PIPO2550,
                     Trees50,
                     PIPO50) %>%
  reduce(full_join, by = "Measurement_ID") %>%
  rbind(none) %>%
  replace(is.na(.), 0) %>% #replace any NA counts with zeros
  group_by(Measurement_ID) %>%
  summarise(Trees2550 = sum(Trees2550, na.rm = T),
            PIPO2550 = sum(PIPO2550, na.rm = T),
            Trees50 = sum(Trees50, na.rm = T),
            PIPO50 = sum(PIPO50, na.rm = T)) %>%
  ungroup()


# Convert to Densities ----------------------------------------------------

tree_data <- tree_tallies %>%
  left_join(meta4, by = "Measurement_ID") %>% #join with dataset with plot size
  mutate(Trees_2550 = Trees2550/md25plot, #convert to densities based on plot sizes
         PIPO_2550 = PIPO2550/md25plot,
         Trees_50 = Trees50/lg50plot,
         PIPO_50 = PIPO50/lg50plot) %>%
  left_join(percPIPO, by = c("Measurement_ID")) %>% #add %pipo variable
  #change NA values of percent pipo to zero when both pipo counts are zero
  mutate(pPIPO = case_when((PIPO_2550 == 0 & PIPO_50 == 0) ~ 0,
                           TRUE ~ pPIPO)) %>%
  dplyr::select(Measurement_ID, #select only measurement ID and density data
                Trees_2550,
                PIPO_2550,
                Trees_50,
                PIPO_50,
                pPIPO)

# Export DF ---------------------------------------------------------------

write.csv(tree_data, here("data_outputs",
                          "01_whwonests",
                               "01_cleaning",
                               "02_nest_tree_densities.csv"))

## END SCRIPT

