
package.list <- c("here", "tidyverse", 
                  "prism", "terra", 'lubridate', 'sf',
                  'st','raster')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Import location dataset -------------------------------------------------
nests <- read.csv(here("data_outputs",
                       "01_whwonests",
                       "01_cleaning",
                       "02_surveys_fates_locations.csv"))

days_nests <- nests %>%
  dplyr::select(Nest_ID, Year_located,
                Visit_date, Prev_date) %>%
  mutate(Visit_date = lubridate::as_date(Visit_date)) %>%
  mutate(Prev_date = lubridate::as_date(Prev_date)) %>%
  mutate(Date = Visit_date) %>%
  group_by(Nest_ID) %>%
  complete(Date = (seq(min(Prev_date),
                             max(Date),
                             by = "1 day"))) %>%
  fill(Year_located, .direction = "updown") %>%
  ungroup() %>%
  dplyr::select(-Visit_date, -Prev_date)

locations <- nests %>%
  distinct(Nest_ID, Year_located, UTM_datum_zone, 
           UTM_N, UTM_E)

# Extraction function -----------------------------------------------------

prism_function <- function(year, variable){
  
  nests2 <- nests %>%
    filter(Year_located == {{year}})
  
  min_date <- nests2 %>%
    filter(Visit_date == min(Visit_date, na.rm = T))%>%
    distinct(Visit_date) %>%
    dplyr::select(Visit_date) %>%
    as_vector()
  
  max_date <- nests2 %>%
    filter(Visit_date == max(Visit_date, na.rm = T))%>%
    dplyr::select(Visit_date) %>%
    distinct(Visit_date) %>%
    as_vector()
  
  get_prism_dailys(type = variable, 
                   minDate = min_date,
                   maxDate = max_date,
                   keepZip = FALSE)
}

#all years for lapply
year_vec <- unique(nests$Year_located)


# Temp --------------------------------------------------------------------

prism_set_dl_dir(path =
                   here('data_raw',
                        '01_whwonests',
                        'climate_data',
                        'prism_daily_tmax'), create = TRUE)

#lapply(year_vec, prism_function, variable = "tmax")

nest_min <- nests %>%
  filter(Prev_date == min(as_date(Prev_date), na.rm = T)) %>%
  distinct(Prev_date) %>%
  dplyr::select(Prev_date) %>%
  as_vector()

nest_max <- nests %>%
  filter(Visit_date == max(Visit_date, na.rm = T)) %>%
  distinct(Visit_date) %>%
  dplyr::select(Visit_date) %>%
  as_vector()

#turn into a raster stack
tmaxstack <- pd_stack(prism_archive_subset("tmax", 
                                         temp_period = "daily",
                                             minDate = nest_min, 
                                             maxDate = nest_max ))

# Prep location data ------------------------------------------------------

#get all the locations ready to combine with climate data
transects10 <- locations %>%
  filter(UTM_datum_zone == "NAD83 10N") %>%
  st_as_sf(coords = c("UTM_E", "UTM_N"), crs = CRS("+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) %>%
  st_transform(crs(raster(tmaxstack)))
transects11 <- locations %>%
  filter(UTM_datum_zone == "NAD83 11N") %>%
  st_as_sf(coords = c("UTM_E", "UTM_N"), crs = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) %>%
  st_transform(crs(raster(tmaxstack)))

transects <- rbind(transects10, transects11); remove(transects10, transects11)

#extract tmax at all transect locations
tmax_locs <- terra::extract(tmaxstack, transects)

#combine transect location with tmax data 
#and get prepped to merge with location data
#so we can filter only dates we care about??
transect_tmax <- transects %>%
  bind_cols(tmax_locs) %>%
  pivot_longer(4:739,
               names_to = "date",
               values_to = "tmax") %>%
  mutate(Climate_Year = str_sub(date, 25, 28)) %>%
  mutate(Climate_Year = as.numeric(Climate_Year)) %>%
  filter(Year_located == Climate_Year) %>%
  mutate(Climate_Month = str_sub(date, 29, 30),
         Climate_Day = str_sub(date, 31, 32),
         Climate_Month = as.numeric(Climate_Month),
         Climate_Day = as.numeric(Climate_Day)) %>%
  unite("Date",
        c(Climate_Year, Climate_Month, Climate_Day),
          sep = "-",
          remove = F) %>%
  mutate(Date = lubridate::as_date(Date))

# PPT ---------------------------------------------------------------------

prism_set_dl_dir(path =
                   here('data_raw',
                        '01_whwonests',
                        'climate_data',
                        'prism_daily_ppt'), create = TRUE)

lapply(year_vec, prism_function, variable = "ppt")


# Extract PPT at each point -----------------------------------------------

#turn into a raster stack
pptstack <- pd_stack(prism_archive_subset("ppt", 
                                           temp_period = "daily",
                                           minDate = nest_min, 
                                           maxDate = nest_max ))

#extract tmax at all transect locations
ppt_locs <- terra::extract(pptstack, transects)

#combine transect location with tmax data 
#and get prepped to merge with location data
#so we can filter only dates we care about??
transect_ppt <- transects %>%
  bind_cols(ppt_locs) %>%
  pivot_longer(4:746,
               names_to = "date",
               values_to = "ppt") %>%
  mutate(Climate_Year = str_sub(date, 24, 27)) %>%
  mutate(Climate_Year = as.numeric(Climate_Year)) %>%
  filter(Year_located == Climate_Year) %>%
  mutate(Climate_Month = str_sub(date, 28, 29),
         Climate_Day = str_sub(date, 30, 31),
         Climate_Month = as.numeric(Climate_Month),
         Climate_Day = as.numeric(Climate_Day)) %>%
  unite("Date",
        c(Climate_Year, Climate_Month, Climate_Day),
        sep = "-",
        remove = F) %>%
  mutate(Date = lubridate::as_date(Date))


# Combine all climate data ------------------------------------------------

days_nests_clim <- days_nests %>%
  left_join(transect_tmax, by = c("Nest_ID",
                                  "Year_located",
                                  "Date")) %>%
  dplyr::select(Nest_ID, Date, Year_located, tmax) %>%
  left_join(transect_ppt, by = c("Nest_ID",
                                 "Year_located",
                                 "Date")) %>%
  dplyr::select(Nest_ID, Date, Year_located, 
                tmax, ppt)

write.csv(days_nests_clim, here('data_outputs',
                              '01_whwonests',
                              '01_cleaning',
                              '03_dailyWeather_nests.csv'))

