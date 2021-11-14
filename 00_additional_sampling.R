#title: Planning for 2021 ARU deployment to fill sampling gaps for CONI
#author: Elly C. Knight
#created: March 17, 2021

library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(rgdal)
library(ggmap)
library(ggspatial)
library(data.table)
library(fasterize)
library(raster)
library(readxl)
library(pgirmess)

#Step 1. Collate sampling locations & metadata----

#ID out missing BU coordinates
aru.bu.missing <- read.csv("Data/APPENDED_REPORT.csv")  %>% 
  mutate(Year=year(ymd(recording_date)),
         SurveyType="ARU") %>% 
  rename(Site=location) %>%
  dplyr::select(SurveyType, Year, Site, longitude, latitude) %>% 
  unique() %>% 
  dplyr::filter(is.na(longitude))
#write.csv(aru.bu.missing, "BU_WLNP_ARUs_missingcoordinates.csv", row.names = FALSE)

#Read everything in & wrangle
aru.wlnp <- read_excel("Data/CONI Database.xlsm", sheet=4) %>% 
  dplyr::filter(!Recognizer=="Boom") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  mutate(Origin="WLNP", SurveyType="ARU") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence)

aru.bu <- read.csv("Data/APPENDED_REPORT.csv") %>% 
  mutate(Year=year(ymd(recording_date)),
         SurveyType="ARU") %>% 
  rename(Site=location) %>% 
  unique() %>% 
  dplyr::filter(!is.na(longitude)) %>% 
  st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>% 
  left_join(read.csv("Data/APPENDED_REPORT.csv") %>% 
              dplyr::filter(species_english_name=="Common Nighthawk") %>% 
              mutate(Presence=1,
                     Year=year(ymd(recording_date))) %>% 
              rename(Site = location) %>% 
              dplyr::select(Year, Site, Presence)) %>% 
  mutate(Presence = ifelse(is.na(Presence), 0, 1),
         Origin="BU") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence) %>% 
  unique()

hum.11 <- read_excel("Data/CONI Database.xlsm", sheet=2) %>% 
  rename(Individuals = '# Individuals') %>% 
  mutate(Year=year(ymd(Date)),
         Presence=ifelse(as.numeric(Individuals) >0, 1, 0)) %>% 
  filter(Zone=="11U",
         !is.na(Presence)) %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26911) %>% 
  st_transform(crs=4326) %>% 
  mutate(Origin="WLNP", SurveyType="Human") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence)

hum.12 <- read_excel("Data/CONI Database.xlsm", sheet=2) %>% 
  rename(Individuals = '# Individuals') %>% 
  mutate(Year=year(ymd(Date)),
         Presence=ifelse(as.numeric(Individuals) >0, 1, 0)) %>% 
  filter(Zone=="12U",
         !is.na(Presence)) %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26912) %>% 
  st_transform(crs=4326) %>% 
  mutate(Origin="WLNP", SurveyType="Human") %>% 
  dplyr::select(Origin, SurveyType, Year, Site, Presence)

#Put together
dat <- rbind(aru.wlnp, aru.bu, hum.11, hum.12) %>% 
  st_coordinates() %>% 
  cbind(rbind(aru.wlnp, aru.bu, hum.11, hum.12)) %>% 
  mutate(surveyID=paste0(Year,"-",SurveyType, "-", Origin))

#Write out
dat.csv <- as.data.frame(dat) %>% 
  dplyr::select(-geometry)
write.csv(dat.csv, "ExistingSamplingLocations.csv", row.names = FALSE)

#1b. Look at repeat sampling----
dat.rep <- dat.csv %>% 
  group_by(X, Y, Site) %>% 
  summarize(n=n()) %>% 
  filter(n > 1) %>% 
  left_join(dat.csv)

#Sites surveyed in 3 years
dat.rep.3 <- dat.rep %>% 
  dplyr::filter(n ==3)

#Sites surveyed at least twice with presence detections
dat.rep.pres <- dat.rep %>% 
  dplyr::filter(Presence==1) %>% 
  group_by(X, Y, Site, n) %>% 
  summarize(Presence=sum(Presence)) %>% 
  mutate(Proportion = Presence/n)

#Step 2. Extract environmental covariates----
dat <- read.csv("ExistingSamplingLocations.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
  st_transform(crs=32612) %>% 
  cbind(read.csv("ExistingSamplingLocations.csv") %>% 
          dplyr::select(X, Y))
 
#2a. Read in and prepare----
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

#Trails
#STILL NEED TO DEAL WITH THIS
trails_official <- read_sf(paste0(gis, "Projects/WLNP/Official_Trails/2015_Official_Trails with AK.shp")) %>% 
  mutate(status="official") %>% 
  dplyr::select(status)
trails_unofficial <- read_sf(paste0(gis, "Projects/WLNP/Unofficial_Trails/UnofficialTrails.shp")) %>% 
  mutate(status="unofficial") %>% 
  dplyr::select(status) 
trails <- rbind(trails_official, trails_unofficial) %>% 
  mutate(status=as.factor(status))
write_sf(trails, "MergedTrails.shp")
trails.sp <- readOGR("MergedTrails.shp")
trails.sp$status <- as.factor(trails.sp$status)
#trails.r <- raster::rasterize(trails.sp, dem.r, field="status")
#writeRaster(rasters/trails.r, "MergedTrails.tiff")

#2b. Extract covs----
coords <- dat %>% 
  st_coordinates()
covs <- data.frame(raster::extract(layers, coords[,c("X", "Y")]))

dat.covs <- cbind(dat, covs) %>% 
  dplyr::filter(!is.na(Soil)) %>% 
  mutate(FireTime=Year-FireHistory,
         FireSeverity = ifelse(Year > 2017, FireSeverity, NA))

write_sf(dat.covs, "ExistingSamplingLocations.shp")

#3. Exploration----

dat.covs <- read_sf("ExistingSamplingLocations.shp")

#Elevation
ggplot(dat.covs) +
  geom_jitter(aes(x=Elevation, y=Presence, colour=factor(Presence))) +
  geom_smooth(aes(x=Elevation, y=Presence))

#TimeSinceFire
ggplot(dat.covs) +
  geom_jitter(aes(x=FireTime, y=Presence, colour=factor(Presence))) +
  geom_smooth(aes(x=FireTime, y=Presence))

#FireSeverity
ggplot(dat.covs) +
  geom_jitter(aes(x=FireSeverity, y=Presence, colour=factor(Presence))) +
  geom_smooth(aes(x=FireSeverity, y=Presence))

#SoilType
soil.dat <- data.frame(table(dat.covs$Presence, dat.covs$Soil)) %>% 
  rename(Presence=Var1, Soil=Var2) %>% 
  pivot_wider(names_from = Presence, names_prefix = "Pres_", values_from="Freq") %>% 
  mutate(Percent = Pres_1/(Pres_0 + Pres_1),
         n=Pres_1+Pres_0)
ggplot(soil.dat) +
  geom_col(aes(x=Soil, y=Percent, fill=n))
levels(soil$PARENT_MAT)

#VegCommunity
veg.dat<- data.frame(table(dat.covs$Presence, dat.covs$Vegetation)) %>% 
  rename(Presence=Var1, Veg=Var2) %>% 
  pivot_wider(names_from = Presence, names_prefix = "Pres_", values_from="Freq") %>% 
  mutate(Percent = Pres_1/(Pres_0 + Pres_1),
         n=Pres_1+Pres_0)
ggplot(veg.dat) +
  geom_col(aes(x=Veg, y=Percent, fill=n))
levels(veg$MAP_BDESC)

#VegetationCover
ggplot(subset(dat.covs, VegetationCover < 4)) +
  geom_jitter(aes(x=VegetationCover, y=Presence, colour=factor(Presence))) +
  geom_smooth(aes(x=VegetationCover, y=Presence))
levels(veg$DENS_MOD)

#VegetationHeight
ggplot(subset(dat.covs, VegetationHeight < 7)) +
  geom_jitter(aes(x=VegetationHeight, y=Presence, colour=factor(Presence))) +
  geom_smooth(aes(x=VegetationHeight, y=Presence))
levels(veg$HT_MOD)

#4. Visualize----
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

map.facet <- ggmap(map_transparent) +
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

map.all <- ggmap(map_transparent) +
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

map.dem <- ggplot() +
  geom_sf(data=trails, aes(colour=status)) +
  geom_spatial_point(aes(x = X, y = Y,
                         colour=factor(Presence),
                         shape=SurveyType),
                     data = dat, 
                     crs=4326,
                     alpha = 0.7,
                     size=2,
                     show.legend = TRUE) +
  scale_colour_viridis_d(name="CONI") +
  facet_grid(SurveyType ~ Year)
map.dem

ggsave(map.dem, file="Figs/Trails_facet.jpeg", width=18, height=6)

#5. Format sampling plan documents----
#Wrangle
locs <- read_sf("/Users/ellyknight/Documents/Employment/Independent Contracts/WatertonCONIMonitoring/GIS/ProposedSamplingLocationsV2.shp") %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  cbind(as.data.frame(read_sf("/Users/ellyknight/Documents/Employment/Independent Contracts/WatertonCONIMonitoring/GIS/ProposedSamplingLocationsV2.shp"))) %>% 
  dplyr::select(-geometry) %>% 
  rename(Group=Group., ObjectiveDetails=Comment, GroupSize=Count) %>% 
  mutate(Objective = ifelse(str_sub(ObjectiveDetails, 1, 8)=="Resample", "Resample", "Gap"),
         ObjectiveDetails = case_when(ObjectiveDetails=="Fire age" ~ "Time since fire",
                                      ObjectiveDetails=="Resample - 0" ~ "Absent before fire",
                                      ObjectiveDetails=="Resample - 1" ~ "Presence before fire",
                                      ObjectiveDetails=="Resample - multiyear" ~ "Multiyear",
                                      ObjectiveDetails=="Soil - 1" ~ "Aeolian sand",
                                      ObjectiveDetails=="Soil - 11" ~ "Excavation",
                                      ObjectiveDetails=="Soil - 13" ~ "Bedrock",
                                      ObjectiveDetails=="Veg - pine" ~ "Lodgepole pine"),
         X=round(X, 6),
         Y=round(Y, 6)) %>% 
  left_join(data.frame(dat) %>%
              dplyr::select(X, Y, Site) %>% 
              mutate(X=round(X, 6),
                     Y=round(Y, 6)) %>% 
              unique()) %>% 
  rename(Latitude=Y, Longitude=X) %>% 
  group_by(Objective) %>% 
  arrange(Latitude) %>% 
  mutate(n=row_number()) %>% 
  ungroup() %>% 
  mutate(Site = ifelse(is.na(Site), paste0("CONI-GAP-", n), Site),
         Group = as.numeric(Group)) %>% 
  dplyr::select(Site, Group, GroupSize, Objective, ObjectiveDetails, Latitude, Longitude)
 
#Prioritize
table(locs$ObjectiveDetails, locs$Group)
#3,4
#2,5
#6,8
#1,7

priorities <- data.frame(Group = c(1:8),
                         Priority = c(8, 3, 1, 2, 4, 5, 7, 6)) %>% 
  arrange(Priority)

locs.prior <- locs %>% 
  left_join(priorities) %>% 
  dplyr::select(-Group) %>% 
  rename(Group=Priority)  %>% 
  dplyr::select(Group, GroupSize, Objective, ObjectiveDetails, Site, Latitude, Longitude) %>% 
  arrange(Group, Site)

#Save as a csv
write.csv(locs.prior, "CommonNightawkSamplingLocations2021.csv", row.names=FALSE)

#Save as a gpx & kml
#CAN'T GET THE SITE NAMES TO WRITE PROPERLY
locs.gpx <- locs.prior %>% 
  dplyr::select(Longitude, Latitude, Site)
coordinates(locs.gpx) <- ~Longitude + Latitude
proj4string(locs.gpx) <- CRS("+proj=longlat +datum=WGS84")

#writeOGR(locs.gpx, dsn="CommonNightawkSamplingLocations2021.gpx", dataset_options="GPX_USE_EXTENSIONS=yes",layer="waypoints",driver="GPX", overwrite_layer = T)

#writeOGR(locs.gpx["Site"], dsn="CommonNightawkSamplingLocations2021.kml", layer="Site", driver="KML", overwrite_layer = T)

#Save as a shp
locs.shp <- locs.prior %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)

write_sf(locs.shp, "CommonNightawkSamplingLocations2021.shp")

#Save out all existing survey locations
locs.all <- dat.covs %>% 
  data.frame() %>% 
  dplyr::select(Site, X, Y) %>% 
  unique()

#Save as csv
write.csv(locs.all, "CommonNightawkSamplingLocationsExisting.csv", row.names=FALSE)

#Save as gpx & kml
coordinates(locs.all) <- ~X+Y
proj4string(locs.all) <- CRS("+proj=longlat +datum=WGS84")

#writeOGR(locs.all, dsn="CommonNightawkSamplingLocationsExisting.gpx", dataset_options="GPX_USE_EXTENSIONS=yes",layer="Site",driver="GPX", overwrite_layer = T)

#writeOGR(locs.all["Site"], dsn="CommonNightawkSamplingLocationsExisting.kml", layer="Site", driver="KML", overwrite_layer = T)

#Save as shp
locs.all.shp <- dat.covs %>% 
  data.frame() %>% 
  dplyr::select(Site, X, Y) %>% 
  unique() %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326)

write_sf(locs.all.shp, "CommonNightawkSamplingLocationsExisting.shp")
