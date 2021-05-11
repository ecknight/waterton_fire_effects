#title: Planning for 2021 ARU deployment to fill sampling gaps for CONI
#author: Elly C. Knight
#created: March 17, 2021

library(tidyverse)
library(lubridate)
library(sf)
library(ggmap)
library(ggspatial)
library(data.table)
library(fasterize)
library(raster)
library(readxl)

#Step 1. Collate sampling locations & metadata----

#ID out missing BU coordinates
aru.bu.missing <- read.csv("PARKS CANADA_Waterton_Lakes_National_Park_2018_2019_2020_report.csv") %>% 
  mutate(Year=year(ymd(recording_date)),
         SurveyType="ARU") %>% 
  rename(Site=location) %>% 
  dplyr::select(SurveyType, Year, Site, longitude, latitude) %>% 
  unique() %>% 
  dplyr::filter(is.na(longitude))
write.csv(aru.bu.missing, "BU_WLNP_ARUs_missingcoordinates.csv", row.names = FALSE)

#Read everything in & wrangle
aru.wlnp <- read.csv("ARUData.csv") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  mutate(Origin="WLNP") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence)

aru.bu <- read.csv("PARKS CANADA_Waterton_Lakes_National_Park_2018_2019_2020_report.csv") %>% 
  mutate(Year=year(ymd(recording_date)),
         SurveyType="ARU") %>% 
  rename(Site=location) %>% 
  unique() %>% 
  dplyr::filter(!is.na(longitude)) %>% 
  st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>% 
  left_join(read.csv("PARKS CANADA_Waterton_Lakes_National_Park_2018_2019_2020_report.csv") %>% 
              dplyr::filter(species_english_name=="Common Nighthawk") %>% 
              mutate(Presence=1,
                     Year=year(ymd(recording_date))) %>% 
              rename(Site = location) %>% 
              dplyr::select(Year, Site, Presence)) %>% 
  mutate(Presence = ifelse(is.na(Presence), 0, 1),
         Origin="BU") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence) %>% 
  unique()

hum.11 <- read.csv("HumanData.csv") %>% 
  filter(Zone=="11U") %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26911) %>% 
  st_transform(crs=4326) %>% 
  mutate(Origin="WLNP") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence)

hum.12 <- read.csv("HumanData.csv") %>% 
  filter(Zone=="12U") %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26912) %>% 
  st_transform(crs=4326) %>% 
  mutate(Origin="WLNP") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence)

#Put together
dat <- rbind(aru.wlnp, aru.bu, hum.11, hum.12) %>% 
  st_coordinates() %>% 
  cbind(rbind(aru.wlnp, aru.bu, hum.11, hum.12)) %>% 
  mutate(surveyID=paste0(Year,"-",SurveyType, "-", Origin))

#Visualize
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

center <- dat  %>% 
  data.frame() %>% 
  summarize(X=mean(X),
            Y=mean(Y))

map <- get_map(center, zoom=11, force=TRUE, maptype="satellite", color="color")

map_attributes <- attributes(map)

map_transparent <- matrix(adjustcolor(map, 
                                      alpha.f = 0.8), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes

dat.map.facet <- ggmap(map_transparent) +
  geom_spatial_point(aes(x = X, y = Y,
                         colour=factor(Presence)),
                     data = dat, 
                     crs=4326,
                     alpha = 0.7,
                     size=2,
                     show.legend = TRUE) +
  facet_grid(SurveyType ~ Year) +
  scale_colour_viridis_d(name="CONI")

ggsave(dat.map.facet, file="OccurrenceVisualization_facet.jpeg", width=14, height=6)

dat.map.all <- ggmap(map_transparent) +
  geom_spatial_point(aes(x = X, y = Y,
                         colour=factor(Presence),
                         shape=SurveyType),
                     data = dat, 
                     crs=4326,
                     alpha = 0.7,
                     size=2,
                     show.legend = TRUE) +
  scale_colour_viridis_d(name="CONI")

ggsave(dat.map.all, file="OccurrenceVisualization_all.jpeg", width=6, height=6)

#Write out
dat.csv <- as.data.frame(dat) %>% 
  dplyr::select(-geometry)
write.csv(dat.csv, "ExistingSamplingLocations.csv", row.names = FALSE)

#Step 2. Extract environmental covariates----
dat <- read.csv("ExistingSamplingLocations.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
  st_transform(crs=32612)
 
#2a. Read in and prepare----
#DEM
dem.r <- raster("/Volumes/ECK004/GIS/Projects/WLNP/DEM_10m.tif") %>% 
  dplyr::filter()
names(dem.r) <- "Elevation"

#Fire history
fire_hist <- read_sf("/Volumes/ECK004/GIS/Projects/WLNP/Fire History/Fire_history2.shp")
fire_hist.r <- fasterize(fire_hist, raster=dem.r, field="YEAR", fun="max")
names(fire_hist.r) <- "FireHistory"

#Kenow fire severity
fire_sev <- read_sf("/Volumes/ECK004/GIS/Projects/WLNP/Kenow 2017 Burn severity/Kenow_severity_classes.shp")
fire_sev.r <- fasterize(fire_sev, raster=dem.r, field="gridcode")
names(fire_sev.r) <- "FireSeverity"

#Soil
soil <- read_sf("/Volumes/ECK004/GIS/Projects/WLNP/Soils_separate_fields/Waterton_Lakes_NP_Soils_1976_Dec_13_2017.shp") %>% 
  mutate(PARENT_MAT=as.factor(PARENT_MAT))
soil.r <- fasterize(soil, raster=dem.r, field="PARENT_MAT")
names(soil.r) <- "Soil"

#Veg
veg <- read_sf("/Volumes/ECK004/GIS/Projects/WLNP/Veg/WATE_VegMap.shp") %>% 
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

#Trails
#STILL NEED TO DEAL WITH THIS
trails_official <- read_sf("/Volumes/ECK004/GIS/Projects/WLNP/Official_Trails/2015_Official_Trails with AK.shp") %>% 
  mutate(status="official") %>% 
  dplyr::select(status)
trails_unofficial <- read_sf("/Volumes/ECK004/GIS/Projects/WLNP/Unofficial_Trails/UnofficialTrails.shp") %>% 
  mutate(status="unofficial") %>% 
  dplyr::select(status) 
trails <- rbind(trails_official, trails_unofficial) %>% 
  mutate(status=as.factor(status))

#2b. Extract covs----
coords <- dat %>% 
  st_coordinates()
covs <- data.frame(raster::extract(layers, coords[,c("X", "Y")]))

dat.covs <- cbind(dat, covs) %>% 
  dplyr::filter(!is.na(Soil))

#3. Exploration----
#Elevation
ggplot(dat.covs) +
  geom_jitter(aes(x=Elevation, y=Presence)) +
  geom_smooth(aes(x=Elevation, y=Presence))

ggplot(dat.covs) +
  geom_jitter(aes(x=FireHistory, y=Presence)) +
  geom_smooth(aes(x=FireHistory, y=Presence))

ggplot(dat.covs) +
  geom_jitter(aes(x=FireSeverity, y=Presence)) +
  geom_smooth(aes(x=FireSeverity, y=Presence))

table(dat.covs$Presence, dat.covs$Soil)
levels(soil$PARENT_MAT)

table(dat.covs$Presence, dat.covs$Vegetation)
levels(veg$MAP_BDESC)

table(dat.covs$Presence, dat.covs$VegetationCover)
levels(veg$DENS_MOD)

table(dat.covs$Presence, dat.covs$VegetationHeight)
levels(veg$HT_MOD)
