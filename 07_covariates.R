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
library(dggridR)

options(scipen=99999)

#1. Wrangle location info----
dat <- read.csv("SurveyDataWithOffsets.csv") %>% 
  mutate(year=year(ymd_hms(DateTime))) %>% 
  group_by(survey, station, year) %>% 
  summarize(detection=max(detection)) %>% 
  ungroup()

#Number of stations
table(dat$year, dat$survey)
table(dat$survey)

#Location data
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

#Put together
dat.locs <- dat %>% 
  mutate(station = gsub(pattern="-00", replacement="-", x=station),
         station = gsub(pattern="-0", replacement="-", x=station),
         station = gsub(pattern="CONI-GAP", replacement="GAP", x=station),
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
 left_join(locs.aru %>% 
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
                     station = gsub(pattern="Highway 6 North 05", replacement="HWY6-5", x=station)))

table(dat.locs$survey)

#Find NAs
dat.na <- dat.locs %>% 
  dplyr::filter(is.na(Longitude))

#6 stations that I can't find utms for

#Create spatial object
dat.sf <- dat.locs %>% 
  dplyr::filter(!station %in% dat.na$station) %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=26912) %>% 
  dplyr::select(survey, station, geometry) %>% 
  unique()

table(dat.sf$survey)

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

#2b. Buffer - trails, soils, veg----
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

#Soil
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
canopy <- veg %>% 
  dplyr::mutate(cover = case_when(DENS_MOD=="1 - Closed Canopy/Continuous (60-100% coverage)" ~ 3,
                                  DENS_MOD=="2 - Open Canopy/Discontinuous (25-60% coverage)" ~ 2,
                                  DENS_MOD=="3 - Dispersed/Sparse Canopy (10-25% coverage)" ~ 1))
summary(canopy$cover)
canopy.r <- fasterize(canopy, raster=dem.r, field="cover")
names(canopy.r) <- "VegetationCover"
writeRaster(canopy.r, "rasters/cover.tif", overwrite=TRUE)

#Veg height
height <- veg %>% 
  dplyr::mutate(height = case_when(HT_MOD=="1 - 30-50 meters (98-162 ft)" ~ 6,
                                   HT_MOD=="2 - 15-30 meters (50-98 feet)" ~ 5,
                                   HT_MOD=="3 - 5-15 meters (16-50 feet)" ~ 4,
                                   HT_MOD=="4 - 0.5-5 meters (1.5-16 feet)" ~ 3,
                                   HT_MOD=="5 - <2 meters" ~ 2,
                                   HT_MOD=="6 - 0.5 meters (< 1.5 feet)" ~ 1))
height.r <- fasterize(height, raster=dem.r, field="height")
names(height.r) <- "VegetationHeight"
writeRaster(height.r, "rasters/height.tif", overwrite=TRUE)

#files <- list.files(path="rasters/", pattern="*.tif")
files <- c("cover.tif", "develop.tif", "grass.tif", "height.tif", "pine.tif", "sand.tif", "trails.tif", "water.tif", "wet.tif", "wetland.tif")
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
    layer.focal <- focal(layer.i, focalWeight(layer.i, d=radius.i, type='circle'), na.rm=TRUE)
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

covs.sf <- dat.sf %>% 
  st_coordinates() %>% 
  cbind(data.frame(dat.sf)) %>% 
  dplyr::select(-geometry) %>% 
  cbind(covs)

table(covs.sf$survey)

write.csv(covs.sf, "SurveyLocationsWithCovs.csv", row.names=FALSE)

#4. Join back to data----
dat <- read.csv("SurveyDataWithOffsets.csv")

dat.covs <- dat %>% 
  mutate(station = gsub(pattern="-00", replacement="-", x=station),
         station = gsub(pattern="-0", replacement="-", x=station),
         station = gsub(pattern="CONI-GAP", replacement="GAP", x=station),
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
  dplyr::filter(!station %in% dat.na$station) %>% 
  left_join(covs.sf) %>% 
  mutate(FireTime=year-FireHistory,
         FireTime=ifelse(FireTime < 0, NA, FireTime),
         FireTime=ifelse(is.na(FireTime), 200, FireTime),
         FireSeverity = ifelse(year > 2017, FireSeverity, NA)) %>% 
  dplyr::filter(!is.na(pine.300),
                Elevation > 0) %>% 
  mutate(ID=paste0(survey,"-",station,"-",year),
         boom=ifelse(detection==2, 1, 0),
         call=ifelse(detection>0, 1, 0),
         DateTime = ymd_hms(DateTime),
         doy = yday(DateTime),
         FireTimeBin = ifelse(FireTime <=4, "y", "n"))

dat.covs.st <- dat.covs %>% 
  dplyr::select(station, survey) %>% 
  unique()

table(dat.covs.st$survey)
#3 less ARU stations that were outside landcover footprint, 1 less that was outside DEM

#5. Thinning----
#Spatial separation of points
dat.dist <- dat.covs %>% 
  dplyr::select(station, X, Y) %>% 
  unique() %>% 
  st_as_sf(coords=c("X", "Y"), crs=26912) %>% 
  st_distance() %>% 
  data.frame() %>% 
  pivot_longer(cols=X1:X268, names_to="station", values_to="distance") %>% 
  mutate(distance = as.integer(distance)) %>% 
  dplyr::filter(distance > 0, distance < 1000)

hist(dat.dist$distance)
#Ok yeah should probably look at spatially thinning

#Set up grid
grid <- dgconstruct(spacing=0.25, metric=TRUE)

dat.station <- dat.covs %>% 
  dplyr::select(station, X, Y) %>% 
  unique() %>% 
  st_as_sf(coords=c("X", "Y"), crs=26912) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(Latitude = Y, Longitude = X) %>% 
  cbind(dat.covs %>% 
          dplyr::select(station, X, Y) %>% 
          unique())

dat.grid <- dat.station %>% 
  mutate(cell = dgGEO_to_SEQNUM(grid, Longitude, Latitude)$seqnum) %>% 
  arrange(cell) %>% 
  full_join(dat.covs)

#Thin
set.seed(1234)
dat.thin <- dat.grid %>% 
  dplyr::select(station, cell, year) %>% 
  unique() %>% 
  group_by(cell, year) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  left_join(dat.grid) 
#pick one survey station per grid cell per year

#Number of sites
dat.grid %>% 
  dplyr::select(ID, survey) %>% 
  unique() %>% 
  group_by(survey) %>% 
  summarize(n = n())

dat.thin %>% 
  dplyr::select(ID, survey) %>% 
  unique() %>% 
  group_by(survey) %>% 
  summarize(n = n())

#could thin more and bootstrap

write.csv(dat.thin, "SurveyDataWithCovs.csv", row.names = FALSE)

#6. VIF----
covs.vif <- read.csv("SurveyDataWithCovs.csv") %>% 
  dplyr::select(Elevation, cover.300, develop.300, grass.300, height.300, pine.300, sand.300, trails.300, water.300, wet.300, wetland.300, FireTime) %>% 
  unique()
M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)
#take out wet

covs.vif <- read.csv("SurveyDataWithCovs.csv") %>% 
  dplyr::select(Elevation, cover.300, develop.300, grass.300, height.300, pine.300, sand.300, trails.300, water.300, wetland.300, FireTime) %>% 
  unique()
M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)
#take out height

covs.vif <- read.csv("SurveyDataWithCovs.csv") %>% 
  dplyr::select(Elevation, cover.300, develop.300, grass.300, pine.300, sand.300, trails.300, water.300, wetland.300, FireTime) %>% 
  unique()

M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)
