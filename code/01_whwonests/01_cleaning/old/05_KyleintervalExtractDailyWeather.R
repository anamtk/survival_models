################################################################################
## This code extracts daily weather data within census intervals from PNW
# transects used to monitor white-headed woodpecker habitat

## Code by Kyle C. Rodman, Ecological Restoration Institute. 
# 4/3/2023

################################################################################
### Bring in necessary packages
package.list <- c("here", "tidyverse", "ncdf4", "raster", "sf")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

################################################################################
### Read in files

## Climate data
climFiles <- list.files(here("Data", "Spatial"), pattern = "_bil.bil", full.names = T)
climFiles <- climFiles[!str_sub(climFiles, start = -7, end = -1) == "aux.xml"]

## Transect locations, reading in as two SF objects (two different UTM zones)
transects10 <- read.csv(here("Data", "Field", "04a_nest_locations_dates.csv"))  %>%
  filter(UTM_datum_zone == "NAD83 10N") %>%
  st_as_sf(coords = c("UTM_E", "UTM_N"), crs = CRS("+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) %>%
  mutate(Julian_start = ifelse(!is.na(Julian_start), Julian_start, Julian_end)) %>% # Infilling NA values for julian start column
  st_transform(crs(raster(climFiles[1])))
transects11 <- read.csv(here("Data", "Field", "04a_nest_locations_dates.csv"))  %>%
  filter(UTM_datum_zone == "NAD83 11N") %>%
  st_as_sf(coords = c("UTM_E", "UTM_N"), crs = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) %>%
  mutate(Julian_start = ifelse(!is.na(Julian_start), Julian_start, Julian_end)) %>% # Infilling NA values for julian start column
  st_transform(crs(raster(climFiles[1])))
transects <- rbind(transects10, transects11); remove(transects10, transects11)

## Check to make sure that all files are in the folder...
sapply(2012:2021, function(i){
  print(length(climFiles[str_detect(climFiles, paste("_",as.character(i), sep =""))]))
})

################################################################################
### Function to get climate data for given location and range of dates

getWeather <- function(pt_index, rastFiles = climFiles, transFiles = transects){
  ## Get the corresponding point in sf object
  pt <- transFiles[pt_index,]
  ## Get only the files from that year - files are consecutive within year
  rastSub <- rastFiles[str_detect(rastFiles, paste("_", pt$Year_located, sep = ""))]
  ## Separate temperature and precipitation
  tmaxGrids <- rastSub[str_detect(rastSub, "tmax")][c(pt$Julian_start:pt$Julian_end)]
  pptGrids <- rastSub[str_detect(rastSub, "ppt")][c(pt$Julian_start:pt$Julian_end)]
  ## Read in as stacks
  tmaxStack <- raster::stack(tmaxGrids)
  pptStack <- raster::stack(pptGrids)
  ## Extract data from these files
  tmax <- raster::extract(tmaxStack, pt); ppt <- raster::extract(pptStack, pt)
  pt$maxTmax_C <- max(tmax, na.rm = T)
  pt$meanTmax_C <- mean(tmax, na.rm = T)
  pt$maxPpt_mm <- max(ppt, na.rm = T)
  pt$meanPpt_mm <- mean(ppt, na.rm = T)
  pt <- st_drop_geometry(pt)
  return(pt)
}

################################################################################
### Apply function to sf object - printing to track progress

## Run through them all
weatherPts <- lapply(1:nrow(transects), function(x){
  ## Track loop
  print(paste("Extracting point", x, "out of", nrow(transects)))
  ## Apply function to row of SF object
  return(getWeather(x))
})

## Merge to one DF
weatherDF <- do.call(rbind, weatherPts)

## Write output to .csv
write_csv(weatherDF, here("Data", "dailyWeather_nests.csv"))

