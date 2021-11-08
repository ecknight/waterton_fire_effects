#title: Covariate extraction for CONI density modelling
#author: Elly C. Knight
#created: November 7, 2021

library(tidyverse)
library(sf)
library(raster)
library(sp)
library(fasterize)
library(lubridate)

options(scipen=99999)

#1. Wrangle location info----
dat <- read.csv("SurveyDataWithOffsets.csv") %>% 
  mutate(year=year(ymd_hms(DateTime))) %>% 
  group_by(survey, station, year) %>% 
  summarize(detection=max(detection)) %>% 
  ungroup()

#Human surveys first
locs.hum <- read.csv("Data/HumanData.csv") %>% 
  rename(station=Site)

dat.hum <- dat %>% 
  dplyr::filter(survey=="human") %>% 
  left_join(locs.hum)

dat.hum.11 <- dat.hum %>% 
  dplyr::filter(Zone=="11U") %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26911) %>% 
  st_transform(crs=4326)
  
dat.hum.12 <- dat.hum %>% 
  dplyr::filter(Zone=="12U") %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26912) %>% 
  st_transform(crs=4326)

dat.hum.sf <- rbind(dat.hum.11, dat.hum.12) %>% 
  dplyr::select(survey, station, year, geometry)

#ARU surveys
locs.2020 <- read.csv("CommonNightawkSamplingLocationsExisting.csv") %>% 
  rename(Latitude=Y, Longitude=X)
locs.2021 <- read.csv("CommonNightawkSamplingLocations2021.csv") %>% 
  dplyr::select(Site, Latitude, Longitude)
locs.bu <- read.csv("Data/BU_WLNP_locations.csv") %>% 
  rename(Site=site, Latitude=Y, Longitude=X) %>% 
  dplyr::select(Site, Longitude, Latitude)
locs.wt <- read.csv("Data/APPENDED_REPORT.csv") %>% 
  rename(Site=location, Latitude=latitude, Longitude=longitude) %>% 
  dplyr::select(Site, Longitude, Latitude)
locs.aru <- rbind(locs.2020, locs.2021, locs.bu, locs.wt) %>% 
  rename(station=Site) %>% 
  unique()

dat.aru.sf <- dat %>% 
  dplyr::filter(survey=="ARU") %>% 
  mutate(station = gsub(pattern="-00", replacement="-", x=station),
         station = gsub(pattern="-0", replacement="-", x=station)) %>% 
  inner_join(locs.aru %>% 
              mutate(station = gsub(pattern="CONI-GAP", replacement="GAP", x=station),
                     station = gsub(pattern="Cardston Gate 0", replacement = "CAR-", x=station),
                     station = gsub(pattern="Cardson Entrance 0", replacement = "ENT-", x=station),
                     station = gsub(pattern="Eskerine 0", replacement = "ESK-", x=station),
                     station = gsub(pattern="HORSESHOE-", replacement = "HORSE-", x=station),
                     station = gsub(pattern="LAKEVIEW-", replacement = "LAKE-", x=station),
                     station = gsub(pattern="OILBASIN", replacement = "OIL", x=station),
                     station = gsub(pattern="Red Rock Parkway ", replacement = "RRP-", x=station),
                     station = gsub(pattern="Belly River 0", replacement = "BEL-", x=station),
                     station = gsub(pattern="Pincher Entrance Transect 0", replacement = "PIN-", x=station),
                     station = gsub(pattern="-00", replacement="-", x=station),
                     station = gsub(pattern="-0", replacement="-", x=station),
                     station = gsub(pattern="Highway 6 North 08", replacement="HWY6-8", x=station),
                     station = gsub(pattern="Highway 6 North 10", replacement="HWY6-10", x=station),
                     station = gsub(pattern="Highway 6 North 11", replacement="HWY6-11", x=station),
                     station = gsub(pattern="Highway 6 South 03", replacement="HWY6-3", x=station),
                     station = gsub(pattern="Highway 6 North 02", replacement="HWY6-2", x=station),
                     station = gsub(pattern="Highway 6 North 05", replacement="HWY6-5", x=station))) %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  dplyr::select(survey, station, year, geometry)

dat.sf <- rbind(dat.aru.sf, dat.hum.sf) %>% 
  st_transform(crs=26912)

#2. Read data and prepare----
gis <- "/Volumes/SSD/GIS/"

#DEM
dem.r <- raster(paste0(gis, "Projects/WLNP/DEM_10m.tif"))
names(dem.r) <- "Elevation"

#Fire history
fire_hist <- read_sf(paste0(gis, "Projects/WLNP/Fire History/Fire_history2.shp"))
fire_hist.r <- fasterize(fire_hist, raster=dem.r, field="YEAR", fun="max")
names(fire_hist.r) <- "FireHistory"

#Kenow fire severity
fire_sev <- read_sf(paste0(gis, "Projects/WLNP/Kenow 2017 Burn severity/Kenow_severity_classes.shp"))
fire_sev.r <- fasterize(fire_sev, raster=dem.r, field="gridcode")
names(fire_sev.r) <- "FireSeverity"

#Soil
soil <- read_sf(paste0(gis, "Projects/WLNP/Soils_separate_fields/Waterton_Lakes_NP_Soils_1976_Dec_13_2017.shp")) %>% 
  mutate(PARENT_MAT=as.factor(PARENT_MAT))
soil.r <- fasterize(soil, raster=dem.r, field="PARENT_MAT")
names(soil.r) <- "Soil"

#Veg
veg <- read_sf(paste0(gis, "/Projects/WLNP/Veg/WATE_VegMap.shp")) %>% 
  mutate(MAP_BDESC=as.factor(MAP_BDESC),
         DENS_MOD=as.factor(DENS_MOD),
         HT_MOD=as.factor(HT_MOD))
veg.r <- fasterize(veg, raster=dem.r, field="MAP_BDESC")
names(veg.r) <- "Vegetation"
canopy.r <- fasterize(veg, raster=dem.r, field="DENS_MOD")
names(canopy.r) <- "VegetationCover"
height.r <- fasterize(veg, raster=dem.r, field="HT_MOD")
names(height.r) <- "VegetationHeight"

layers <- stack(dem.r, fire_hist.r, fire_sev.r, soil.r, veg.r, canopy.r, height.r)

#3. Extract----
coords <- dat.sf %>% 
  st_coordinates()
covs <- data.frame(raster::extract(layers, coords[,c("X", "Y")]))

dat.covs <- dat.sf %>% 
  st_coordinates() %>% 
  cbind(data.frame(dat.sf)) %>% 
  dplyr::select(-geometry) %>% 
  cbind(covs) %>% 
  dplyr::filter(!is.na(Soil)) %>% 
  mutate(FireTime=year-FireHistory,
         FireTime=ifelse(FireTime < 0, NA, FireTime),
         FireSeverity = ifelse(year > 2017, FireSeverity, NA)) 

write.csv(dat.covs, "SurveyLocationsWithCovs.csv", row.names=FALSE)
