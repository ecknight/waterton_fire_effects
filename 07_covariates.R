#title: Covariate extraction for CONI density modelling
#author: Elly C. Knight
#created: November 7, 2021

library(tidyverse)
library(sf)
library(raster)
library(sp)
library(fasterize)
library(lubridate)
library(corrplot)
library(usdm)

options(scipen=99999)

#1. Wrangle location info----
dat <- read.csv("SurveyDataWithOffsets.csv") %>% 
  mutate(year=year(ymd_hms(DateTime))) %>% 
  group_by(survey, station, year) %>% 
  summarize(detection=max(detection)) %>% 
  ungroup()

#ARU surveys
locs.2020 <- read.csv("CommonNightawkSamplingLocationsExisting.csv") %>% 
  rename(Latitude=Y, Longitude=X) %>% 
  group_by(Site) %>% 
  sample_n(1) %>% #Randomly sample 1 of each of the highway 6 south 4-6 stations
  ungroup()

locs.2021 <- read.csv("CommonNightawkSamplingLocations2021.csv") %>% 
  dplyr::filter(Objective=="Gap") %>%  #only include new stations from 2021
  dplyr::select(Site, Latitude, Longitude)

locs.wt <- read.csv("Data/PARKS CANADA_locations_202119_WLNPfilter.csv") %>% 
  rename(Site = location, Latitude=latitude, Longitude=longitude) %>% 
  dplyr::filter(!Site %in% locs.2020$Site,
                !Site=="WLNP-SUM-09") %>% #remove old WLNP-SUM-09 coordinates
  dplyr::select(Site, Longitude, Latitude)

locs.aru <- rbind(locs.2020, locs.2021, locs.wt) %>% 
  rename(station=Site) %>% 
  mutate(Longitude = round(as.numeric(Longitude), digits=5),
         Latitude = round(as.numeric(Latitude), digits=5)) %>% 
  unique()

dat.sf <- dat %>% 
  mutate(station = gsub(pattern="-00", replacement="-", x=station),
         station = gsub(pattern="-0", replacement="-", x=station)) %>% 
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
         station = gsub(pattern="Highway 6 North 05", replacement="HWY6-5", x=station)) %>% 
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
  st_transform(crs=26912) %>% 
  dplyr::select(station, geometry) %>% 
  unique()

#2. Read data and prepare----
gis <- "/Volumes/SSD/GIS/"

#2a. Point level - DEM, fire----
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

#2b. Buffer - trails, soils, veg
#Trails
trails.r <- raster("rasters/MergedTrails.tif")
plot(trails.r)
rclmat <- matrix(c(1, 2, 1, 
                   NA, NA, 0), 
                 ncol=3, byrow=TRUE) 
trails <- reclassify(trails.r, rclmat)
plot(trails)
names(trails) <- "trails"
writeRaster(trails, "rasters/trails.tif", overwrite=TRUE)

#Soil - NOT SURE WHAT TO DO WITH THIS####
soil <- read_sf(paste0(gis, "Projects/WLNP/Soils_separate_fields/Waterton_Lakes_NP_Soils_1976_Dec_13_2017.shp")) %>% 
  mutate(PARENT_MAT=as.factor(PARENT_MAT))
soil.r <- fasterize(soil, raster=dem.r, field="PARENT_MAT")
names(soil.r) <- "Soil"
levels(soil$PARENT_MAT)
rclmat <- matrix(c(0, 2.5, 1, 
                   2.5, 5.5, 0,
                   5.5, 6.5, 0,
                   6.5, 11.5, 0,
                   11.5, 12.5, 0,
                   12.5, 13.5, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
sand <- reclassify(soil.r, rclmat)
plot(sand)
names(sand) <- "sand"
writeRaster(sand, "rasters/sand.tif", overwrite=TRUE)

#Veg
veg <- read_sf(paste0(gis, "/Projects/WLNP/Veg/WATE_VegMap.shp")) %>% 
  mutate(MAP_BDESC=as.factor(MAP_BDESC),
         DENS_MOD=as.factor(DENS_MOD),
         HT_MOD=as.factor(HT_MOD))
veg.r <- fasterize(veg, raster=dem.r, field="MAP_BDESC")
names(veg.r) <- "Vegetation"
levels(veg$MAP_BDESC)

#Water
rclmat <- matrix(c(0, 24.5, 0, 
                   24.5, 25.5, 1,
                   25.5, 30.5, 0,
                   30.5, 31.5, 1,
                   31.5, 37.5, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
water <- reclassify(veg.r, rclmat)
plot(water)
names(water) <- "water"
writeRaster(water, "rasters/water.tif", overwrite=TRUE)

#Wetland
rclmat <- matrix(c(0, 12.5, 0, 
                   12.5, 14.5, 1,
                   14.5, 25.5, 0,
                   25.5, 26.5, 1,
                   26.5, 29.5, 0,
                   29.5, 30.5, 1,
                   30.5, 34.5, 0,
                   34.5, 35.5, 1,
                   35.5, 37.5, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
wetland <- reclassify(veg.r, rclmat)
plot(wetland)
names(wetland) <- "wetland"
writeRaster(wetland, "rasters/wetland.tif", overwrite=TRUE)

#Wet
rclmat <- matrix(c(0, 24.5, 0, 
                   24.5, 26.5, 1,
                   26.5, 29.5, 0,
                   29.5, 31.5, 1,
                   31.5, 34.5, 0,
                   34.5, 35.5, 1,
                   35.5, 37.5, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
wet <- reclassify(veg.r, rclmat)
plot(wet)
names(wet) <- "wet"
writeRaster(wet, "rasters/wet.tif", overwrite=TRUE)

#Pine
rclmat <- matrix(c(0, 17.5, 0, 
                   17.5, 20.5, 1,
                   20.5, 37.5, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
pine <- reclassify(veg.r, rclmat)
plot(pine)
names(pine) <- "pine"
writeRaster(pine, "rasters/pine.tif", overwrite=TRUE)

#Grassland
rclmat <- matrix(c(0, 15.5, 0, 
                   15.5, 16.5, 1,
                   16.5, 37.5, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
grass <- reclassify(veg.r, rclmat)
plot(grass)
names(grass) <- "grass"
writeRaster(grass, "rasters/grass.tif", overwrite=TRUE)

#Developed
rclmat <- matrix(c(0, 27.5, 0, 
                   27.5, 29.5, 1,
                   29.5, 37.5, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
develop <- reclassify(veg.r, rclmat)
plot(develop)
names(develop) <- "develop"
writeRaster(develop, "rasters/develop.tif", overwrite=TRUE)

#Veg cover
canopy.r <- fasterize(veg, raster=dem.r, field="DENS_MOD")
names(canopy.r) <- "VegetationCover"
writeRaster(canopy.r, "rasters/cover.tif", overwrite=TRUE)

#Veg height
height.r <- fasterize(veg, raster=dem.r, field="HT_MOD")
names(height.r) <- "VegetationHeight"
writeRaster(height.r, "rasters/height.tif", overwrite=TRUE)

files <- list.files(path="rasters/", pattern="*.tif")
files <- "sand.tif"
radii <- c(300)
loop <- expand.grid(files=files, radius=radii)

for(i in 1:nrow(loop)){
  name.i <- str_sub(loop$files[i], -100, -5)
  layer.i <- raster(paste0("rasters/", as.character(loop$files[i])))
  radius.i <- loop$radius[i]
  if(name.i %in% c("cover", "height")){
    layer.focal <- focal(layer.i, focalWeight(layer.i, d=radius.i, type='circle'), na.rm=TRUE)
  }
  else{
    layer.focal <- focal(layer.i, focalWeight(layer.i, d=radius.i, type='circle'))
  }
  names(layer.focal) <- paste0(name.i,"-",radius.i)
  writeRaster(layer.focal, paste0("rasters/",name.i,"-",radius.i,".tif"), overwrite=TRUE)
  print(paste0("Completed raster ", name.i, "-", radius.i, " - ", i, " of ", nrow(loop), " rasters"))
}

cover.300 <- raster("rasters/cover-300.tif")
develop.300 <- raster("rasters/develop-300.tif")
grass.300 <- raster("rasters/grass-300.tif")
height.300 <- raster("rasters/height-300.tif")
trails.300 <- raster("rasters/trails-300.tif")
pine.300 <- raster("rasters/pine-300.tif")
sand.300 <- raster("rasters/sand-300.tif")
water.300 <- raster("rasters/water-300.tif")
wet.300 <- raster("rasters/wet-300.tif")
wetland.300 <- raster("rasters/wetland-300.tif")

layers <- stack(dem.r, fire_hist.r, fire_sev.r, cover.300, develop.300, grass.300, height.300, trails.300, pine.300, sand.300, water.300, wet.300, wetland.300)

#3. Extract----
coords <- dat.sf %>% 
  st_coordinates()
covs <- data.frame(raster::extract(layers, coords[,c("X", "Y")]))

dat.covs <- dat.sf %>% 
  st_coordinates() %>% 
  cbind(data.frame(dat.sf)) %>% 
  dplyr::select(-geometry) %>% 
  cbind(covs)

write.csv(dat.covs, "SurveyLocationsWithCovs.csv", row.names=FALSE)

#4. VIF----
covs.vif <- dat.covs %>% 
  dplyr::select(Elevation, cover.300, develop.300, grass.300, height.300, trails.300, pine.300, sand.300, water.300, wet.300, wetland.300)
M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)
#take out wet

covs.vif <- dat.covs %>% 
  dplyr::select(Elevation, cover.300, develop.300, grass.300, height.300, trails.300, pine.300, sand.300, water.300, wetland.300)
M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)
#take out height

covs.vif <- dat.covs %>% 
  dplyr::select(Elevation, cover.300, develop.300, grass.300, trails.300, pine.300, sand.300, water.300, wetland.300)
M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)