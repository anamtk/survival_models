# Prep landscape variables for nests
# Ana Miller-ter Kuile
# June 6, 2022

# this script preps the landscape variables for nests, 
# only at the 1000m radius scale. Thus it's a little
# less scary than the other script 02c_...

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


# Load data ---------------------------------------------------------------

emfwor <- read.csv(here("data_raw",
                        "landscape_data",
                        "EMFWOR_LF-based_IPM_NestSurv_data_n68_nests.csv"))

emmaor <- read.csv(here("data_raw",
                        "landscape_data", 
                        "EMMAOR_LF-based_IPM_NestSurv_data_n149_nests.csv"))

empaid <- read.csv(here("data_raw",
                        "landscape_data", 
                        "EMPAID_LF-based_IPM_NestSurv_data_n146_nests.csv"))

tx <- read.csv(here('data_raw',
                    'landscape_data',
                    "IPM_Nests_BurnTypeSeasonMar2023.csv"))

land <- emfwor %>%
  rbind(emmaor) %>%
  rbind(empaid)
# Select variables --------------------------------------------------------

colnames(land)
 
#taking just closed canopy forest variables from that set of variables,
#along with treatment variables
forest <- land %>%
  dplyr::select(FID, Nest_ID, Year_located, Date_located, Project_ID,
                Transect_ID, Point_ID, Pt_Trt,
                X3f1areaam1, X3f1areaam2,      
                X3f1areacv1, X3f1areacv2,
                X3f1areamn1, X3f1areamn2,      
                X3f1ca1, X3f1ca2,           
                X3f1contag, X3f1lpi, X3f1lpi1,          
                X3f1lpi2, X3f1np, X3f1np1,           
                X3f1np2, X3f1pd, X3f1pd1,           
                X3f1pd2,  X3f1pland1,  X3f1pland2,        
                X3f1proxam, X3f1proxam1, X3f1proxam2,       
                X3f1proxcv, X3f1proxcv1, X3f1proxcv2,         
                X3f1proxmn, X3f1proxmn1, X3f1proxmn2,
                X9f1areaam1, X9f1areaam2,       
                X9f1areacv1, X9f1areacv2, 
                X9f1areamn1, X9f1areamn2,       
                X9f1ca1, X9f1ca2 ,          
                X9f1contag, X9f1lpi, X9f1lpi1,          
                X9f1lpi2, X9f1np, X9f1np1,           
                X9f1np2, X9f1pd, X9f1pd1,           
                X9f1pd2, X9f1pland1, X9f1pland2,        
                X9f1proxam, X9f1proxam1, X9f1proxam2,       
                X9f1proxcv, X9f1proxcv1, X9f1proxcv2,       
                X9f1proxmn, X9f1proxmn1, X9f1proxmn2,
                Ha11_1000, Ha12_1000, Ha13_1000, Ha14_1000,
                Ha15_1000, Ha16_1000, Ha17_1000, Ha18_1000,
                Ha19_1000, Ha20_1000, Ha21_1000) %>%
  #rename with spatial scale and year for 2013 and 2019 data
  rename_with(~str_replace(., pattern = "X3f1", replacement = "a1000_2013_"), 
              contains('X3f1')) %>%
  rename_with(~str_replace(., pattern = "X9f1", replacement = "a1000_2019_"),
              contains("X9f1"))
       
#get burn treatment new categories       
tx2 <- tx %>%
  dplyr::select(Nest_ID, 
                P11_1000, P12_1000, P13_1000, P14_1000,
                P15_1000, P16_1000, P17_1000, P18_1000,
                P19_1000, P20_1000, P21_1000,
                Rx11_1000, Rx12_1000, Rx13_1000, Rx14_1000,
                Rx15_1000, Rx16_1000, Rx17_1000, Rx18_1000,
                Rx19_1000, Rx20_1000, Rx21_1000,
                WF11_1000, WF12_1000, WF13_1000, WF14_1000,
                WF15_1000, WF16_1000, WF17_1000, WF18_1000,
                WF19_1000, WF20_1000, WF21_1000)
                
#combine these two datasets
forest2 <- forest %>%
  left_join(tx2, by = 'Nest_ID')

forest3 <- forest2 %>%
  mutate_if(., is.integer, as.numeric) %>%
  #weighted mean forest patch size
  mutate(a1000_areaam1 = case_when(Year_located <= 2013 ~ a1000_2013_areaam1,
                                 Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_areaam1,
                                 Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_areaam1,
                                 TRUE ~ NA_real_)) %>%
  mutate(a1000_areaam2 = case_when(Year_located <= 2013 ~ a1000_2013_areaam2,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_areaam2,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_areaam2,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_areacv1 = case_when(Year_located <= 2013 ~ a1000_2013_areacv1,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_areacv1,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_areacv1,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_areacv2 = case_when(Year_located <= 2013 ~ a1000_2013_areacv2,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_areacv2,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_areacv2,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_areamn1 = case_when(Year_located <= 2013 ~ a1000_2013_areamn1,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_areamn1,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_areamn1,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_areamn2 = case_when(Year_located <= 2013 ~ a1000_2013_areamn2,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_areamn2,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_areamn2,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_ca1 = case_when(Year_located <= 2013 ~ a1000_2013_ca1,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_ca1,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_ca1,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_ca2 = case_when(Year_located <= 2013 ~ a1000_2013_ca2,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_ca2,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_ca2,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_contag = case_when(Year_located <= 2013 ~ a1000_2013_contag,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_contag,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_contag,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_lpi = case_when(Year_located <= 2013 ~ a1000_2013_lpi,
                                   Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_lpi,
                                   Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_lpi,
                                   TRUE ~ NA_real_)) %>%
  mutate(a1000_lpi1 = case_when(Year_located <= 2013 ~ a1000_2013_lpi1,
                               Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_lpi1,
                               Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_lpi1,
                               TRUE ~ NA_real_)) %>%
  mutate(a1000_lpi2 = case_when(Year_located <= 2013 ~ a1000_2013_lpi2,
                               Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_lpi2,
                               Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_lpi2,
                               TRUE ~ NA_real_)) %>%
  mutate(a1000_np = case_when(Year_located <= 2013 ~ a1000_2013_np,
                               Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_np,
                               Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_np,
                               TRUE ~ NA_real_)) %>%
  mutate(a1000_np1 = case_when(Year_located <= 2013 ~ a1000_2013_np1,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_np1,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_np1,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_np2 = case_when(Year_located <= 2013 ~ a1000_2013_np2,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_np2,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_np2,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_pd = case_when(Year_located <= 2013 ~ a1000_2013_pd,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_pd,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_pd,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_pd1 = case_when(Year_located <= 2013 ~ a1000_2013_pd1,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_pd1,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_pd1,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_pd2 = case_when(Year_located <= 2013 ~ a1000_2013_pd2,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_pd2,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_pd2,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_pland1 = case_when(Year_located <= 2013 ~ a1000_2013_pland1,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_pland1,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_pland1,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_pland2 = case_when(Year_located <= 2013 ~ a1000_2013_pland2,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_pland2,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_pland2,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_proxam = case_when(Year_located <= 2013 ~ a1000_2013_proxam,
                              Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxam,
                              Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxam,
                              TRUE ~ NA_real_)) %>%
  mutate(a1000_proxam1 = case_when(Year_located <= 2013 ~ a1000_2013_proxam1,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxam1,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxam1,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_proxam2 = case_when(Year_located <= 2013 ~ a1000_2013_proxam2,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxam2,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxam2,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_proxcv = case_when(Year_located <= 2013 ~ a1000_2013_proxcv,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxcv,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxcv,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_proxcv1 = case_when(Year_located <= 2013 ~ a1000_2013_proxcv1,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxcv1,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxcv1,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_proxcv2 = case_when(Year_located <= 2013 ~ a1000_2013_proxcv2,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxcv2,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxcv2,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_proxmn = case_when(Year_located <= 2013 ~ a1000_2013_proxmn,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxmn,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxmn,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_proxmn1 = case_when(Year_located <= 2013 ~ a1000_2013_proxmn1,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxmn1,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxmn1,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_proxmn2 = case_when(Year_located <= 2013 ~ a1000_2013_proxmn2,
                                  Year_located > 2013 & Pt_Trt != "U" ~ a1000_2019_proxmn2,
                                  Year_located > 2013 & Pt_Trt == "U" ~ a1000_2013_proxmn2,
                                  TRUE ~ NA_real_)) %>%
  mutate(a1000_Ha = case_when(Year_located == 2012 ~ Ha12_1000,
                              Year_located == 2013 ~ Ha13_1000,
                              Year_located == 2014 ~ Ha14_1000,
                              Year_located == 2015 ~ Ha15_1000,
                              Year_located == 2016 ~ Ha16_1000,
                              Year_located == 2017 ~ Ha17_1000,
                              Year_located == 2018 ~ Ha18_1000,
                              Year_located == 2019 ~ Ha19_1000,
                              Year_located == 2020 ~ Ha20_1000,
                              Year_located == 2021 ~ Ha21_1000,
                              TRUE ~ NA_real_)) %>%
  #pile burned area (probs not using)
  mutate(a1000_PBu = case_when(Year_located == 2012 ~ P12_1000,
                              Year_located == 2013 ~ P13_1000,
                              Year_located == 2014 ~ P14_1000,
                              Year_located == 2015 ~ P15_1000,
                              Year_located == 2016 ~ P16_1000,
                              Year_located == 2017 ~ P17_1000,
                              Year_located == 2018 ~ P18_1000,
                              Year_located == 2019 ~ P19_1000,
                              Year_located == 2020 ~ P20_1000,
                              Year_located == 2021 ~ P21_1000,
                              TRUE ~ NA_real_)) %>%
  #burned area
  mutate(a1000_RxBu = case_when(Year_located == 2012 ~ Rx12_1000,
                                Year_located == 2013 ~ Rx13_1000,
                                Year_located == 2014 ~ Rx14_1000,
                                Year_located == 2015 ~ Rx15_1000,
                                Year_located == 2016 ~ Rx16_1000,
                                Year_located == 2017 ~ Rx17_1000,
                                Year_located == 2018 ~ Rx18_1000,
                                Year_located == 2019 ~ Rx19_1000,
                                Year_located == 2020 ~ Rx20_1000,
                                Year_located == 2021 ~ Rx21_1000,
                                TRUE ~ NA_real_)) %>%
  #wildfire burned area
  mutate(a1000_WFBu = case_when(Year_located == 2012 ~ WF12_1000,
                               Year_located == 2013 ~ WF13_1000,
                               Year_located == 2014 ~ WF14_1000,
                               Year_located == 2015 ~ WF15_1000,
                               Year_located == 2016 ~ WF16_1000,
                               Year_located == 2017 ~ WF17_1000,
                               Year_located == 2018 ~ WF18_1000,
                               Year_located == 2019 ~ WF19_1000,
                               Year_located == 2020 ~ WF20_1000,
                               Year_located == 2021 ~ WF21_1000,
                               TRUE ~ NA_real_)) %>%
  #select all the final variables
  dplyr::select(Nest_ID, 
                a1000_areaam1, a1000_areaam2,
                a1000_areacv1, a1000_areacv2,
                a1000_areamn1, a1000_areamn2,
                a1000_ca1, a1000_ca2,
                a1000_contag, 
                a1000_lpi,a1000_lpi1, a1000_lpi2,
                a1000_np, a1000_np1, a1000_np2,
                a1000_pd, a1000_pd1, a1000_pd2,
                a1000_pland1, a1000_pland2,
                a1000_proxam, a1000_proxam1, a1000_proxam2,
                a1000_proxcv, a1000_proxcv1, a1000_proxcv2,
                a1000_proxmn, a1000_proxmn1, a1000_proxmn2,
                a1000_Ha ,a1000_PBu, a1000_RxBu, a1000_WFBu)

# Export ------------------------------------------------------------------

write.csv(forest3, here("data_outputs",
                        "01_cleaning",
                        "03_nest_landscape_1000_covariates.csv"))

#END SCRIPT
