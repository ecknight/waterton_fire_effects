library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(raster)
library(ggmap)
library(ggspatial)
library(patchwork)
library(gridExtra)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))


#FIGURE 1. STUDY AREA----
dat <- read.csv("SurveyDataWithOffsets&Covariates.csv") %>% 
  mutate(Kenow = ifelse(FireHistory==2017, "impact", "control"),
         Kenow = ifelse(is.na(Kenow), "control", Kenow)) %>% 
  dplyr::select(survey, station, year, X, Y) %>% 
  unique()

#1a. Study area----
nam <- map_data("world", region=c("Canada", 
                                  "USA", 
                                  "Mexico",
                                  "Guatemala", 
                                  "Belize", 
                                  "El Salvador",
                                  "Honduras", 
                                  "Nicaragua", 
                                  "Costa Rica",
                                  "Panama", 
                                  "Jamaica", 
                                  "Cuba", 
                                  "The Bahamas",
                                  "Haiti", 
                                  "Dominican Republic", 
                                  "Antigua and Barbuda",
                                  "Dominica", 
                                  "Saint Lucia", 
                                  "Saint Vincent and the Grenadines", 
                                  "Barbados",
                                  "Grenada",
                                  "Trinidad and Tobago")) %>% 
  dplyr::filter(!group%in%c(258:264))

nam.eq <- nam %>% 
  st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>%
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(nam)

area.eq.center <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=26912) %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  dplyr::summarize(X=mean(X),
                   Y=mean(Y))

map.nam <- ggplot() +
  geom_polygon(data=nam.eq, aes(x=X, y=Y, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_text(label="★", aes(x=X, y=Y), size=5, family = "HiraKakuPro-W3", data=area.eq.center, colour="black") +
  xlim(c(-4000000, 3000000)) +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  my.theme +
  theme(plot.margin = unit(c(0,0,-1,-1), "cm"))
map.nam

#Part 2. Study sites----

#Locations
dat.wgs <-dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=26912) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  rename(lon=X, lat=Y) %>% 
  cbind(dat)

area.wgs.center <- dat.wgs %>% 
  dplyr::summarize(X=mean(lon),
                   Y=mean(lat))
area.wgs.center


#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

map <- get_map(area.wgs.center, zoom=10, force=TRUE, maptype="satellite", color="color")

map_attributes <- attributes(map)

map_transparent <- matrix(adjustcolor(map, 
                                      alpha.f = 0.7), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes

#Get park perimeter
wlnp <- read_sf("/Volumes/SSD/GIS/Administrative/Canada/National Parks/CLAB_AB_2021-12-02/CLAB_AB_2021-12-02.shp") %>% 
  dplyr::filter(CLAB_ID=="WATE") %>% 
  st_transform(crs=4326)

wlnp.sp <- as_Spatial(wlnp)

#Get fire perimeter
kenow <- read_sf("/Volumes/SSD/GIS/Projects/WLNP/Kenow 2017 Burn severity/Kenow_severity_classes.shp") %>% 
#  st_union() %>% 
  st_transform(crs=4326)
kenow.sp <- as_Spatial(kenow)

fire <- read_sf("/Volumes/SSD/GIS/Projects/WLNP/Fire History/Fire_history2.shp") %>% 
  dplyr::filter(YEAR==2017) %>% 
  st_transform(crs=4326) %>% 
  st_make_valid() %>% 
  st_union()

fire.sp <- as_Spatial(fire)

#map
map.site <- ggmap(map_transparent) +
  geom_polygon(data=wlnp.sp, aes(x=long, y=lat), fill=NA, colour="black") +
#  geom_polygon(data=kenow.sp, aes(x=long, y=lat), fill=NA, colour="red", alpha = 0.4) +
#  geom_polygon(data=fire.sp, aes(x=long, y=lat), fill="red", alpha=0.4) +
  geom_spatial_point(aes(x = lon, y = lat,
                         colour=factor(year),
                         shape=survey),
                     data = dat.wgs, 
                     crs=4326,
                     alpha = 0.8,
                     size=4) +
  scale_shape_manual(values=c(16, 17), name="Survey type") +
  scale_colour_viridis_d(name="Survey year") +
  ggspatial::annotation_north_arrow(location = "tl",
                                    style = ggspatial::north_arrow_orienteering(fill = c("grey80", "grey20"), line_col = "grey20")) +
  ggsn::scalebar(x.min = -114.2, x.max = -114.05, 
                 y.min = 48.96, y.max = 49.05, 
                 transform=TRUE, model="WGS84",
                 dist=5, dist_unit="km",
                 box.fill=c("grey80", "grey20"),
                 box.color="grey20",
                 height = 0.15, 
                 st.dist = 0.1) +
  coord_sf(crs=4326) +
  my.theme +
  xlab("") +
  ylab("") +
  xlim(c(-114.2, -113.6)) +
  ylim(c(48.95, 49.25)) +
  theme(legend.position = "right",
        axis.text.x.bottom=element_text(size=10),
        axis.text.y.left=element_text(size=10))
#map.site

#1c. Put it together####
plot.sa <- map.site +
  inset_element(map.nam,
                right=0.98,
                bottom=0.55,
                left=0.70,
                top=0.98)
#plot.sa

ggsave(plot.sa, filename="Figs/StudyArea.jpeg", device="jpeg", width=9, height=6, units="in", dpi=600)

#FIGURE 4. OFFSETS-----
pr <- read.csv("OffsetPredictionsForFigure.csv") %>%
  mutate(Date=as.Date(jday*365, origin="2015-01-01"))

pr.boom <- pr %>% 
  dplyr::filter(!is.na(p.boom))

dat.boom <- read.csv("OffsetDataForFigure_Boom.csv") %>% 
  mutate(doy = yday(ymd_hms(DateTime)),
         Date = as.Date(doy, origin="2015-01-01")) %>% 
  dplyr::filter(doy >= min(pr.boom$doy),
                doy <= max(pr.boom$doy),
                tsss >= min(pr.boom$tsss),
                tsss <= max(pr.boom$tsss))

dat.call <- read.csv("OffsetDataForFigure_Call.csv") %>% 
  mutate(doy = yday(ymd_hms(DateTime)),
         Date = as.Date(doy, origin="2015-01-01")) %>% 
  dplyr::filter(doy >= min(pr$doy),
                doy <= max(pr$doy),
                tsss >= min(pr$tsss),
                tsss <= max(pr$tsss))

dat <- read.csv("SurveyDataWithOffsets.csv") %>% 
  mutate(doy = yday(ymd_hms(DateTime)),
         Date = as.Date(doy, origin="2015-01-01")) %>% 
  dplyr::filter(survey=="ARU")

dat.boom <- dat %>% 
  dplyr::filter(doy >= min(pr.boom$doy),
                doy <= max(pr.boom$doy),
                tsss >= min(pr.boom$tsss),
                tsss <= max(pr.boom$tsss)) %>% 
  mutate(detection = ifelse(detection==2, 1, 0))

dat.call <- dat %>% 
  dplyr::filter(doy >= min(pr$doy),
                doy <= max(pr$doy),
                tsss >= min(pr$tsss),
                tsss <= max(pr$tsss)) %>% 
  mutate(detection = ifelse(detection %in% c(1,2), 1, 0))

plot.boom.offset <- ggplot() +
  geom_raster(aes(x=Date, y=tsss, fill=p.boom), data=pr.boom, alpha=0.7) +
  scale_fill_viridis_c(name="Probability\nof territorial\nactivity", direction=-1) +
  geom_point(aes(x=Date, y=tsss, colour=factor(detection)), alpha=0.5, data=subset(dat.boom, detection==0)) +
  geom_point(aes(x=Date, y=tsss, colour=factor(detection)), data=subset(dat.boom, detection==1)) +
  scale_colour_manual(values=c("grey70", "grey30"), name="Common\nnighthawk\ndetection", labels=c("Absent", "Present")) +
  ylab("Hours since sunset") +
  my.theme
plot.boom.offset

plot.call.offset <- ggplot() +
  geom_raster(aes(x=Date, y=tsss, fill=p.call), data=pr, alpha=0.7) +
  scale_fill_viridis_c(name="Probability\nof home range\nactivity", direction=-1) +
  geom_point(aes(x=Date, y=tsss, colour=factor(detection)), alpha=0.5, data=subset(dat.call, detection==0)) +
  geom_point(aes(x=Date, y=tsss, colour=factor(detection)), data=subset(dat.call, detection==1)) +
  scale_colour_manual(values=c("grey70", "grey30"), name="Common\nnighthawk\ndetection", labels=c("Absent", "Present")) +
  ylab("Hours since sunset") +
  my.theme
plot.call.offset

ggsave(plot=grid.arrange(plot.boom.offset, plot.call.offset, ncol=2), filename="Figs/Figure4Offsets.jpeg", device="jpeg", width=12, height=5, units="in")

#FIGURE 5. COVARIATE EFFECTS####
newdat <- read.csv("OccupancyModelPredictions.csv")

newdat.grass <- newdat %>% 
  dplyr::filter(elevation==1570)

newdat.elev <- newdat %>% 
  dplyr::filter(grass==0.14)

plot.boom.grass <- ggplot(newdat.grass) +
  geom_ribbon(aes(x=grass, ymin=delta.boom.low, ymax = delta.boom.high), alpha = 0.5) +
  geom_line(aes(x=grass, y=delta.boom), show.legend = FALSE) +
  xlab("") +
  ylab("Territorial habitat suitability") +
  my.theme

plot.call.grass <- ggplot(newdat.grass) +
  geom_ribbon(aes(x=grass, ymin=delta.call.low, ymax = delta.call.high), alpha = 0.5) +
  geom_line(aes(x=grass, y=delta.call)) +
  xlab("Proportion of grassland in 300 m") +
  ylab("Home range habitat suitability") +
  my.theme

plot.boom.elev <- ggplot(newdat.elev) +
  geom_ribbon(aes(x=elevation, ymin=delta.boom.low, ymax = delta.boom.high), alpha = 0.5) +
  geom_line(aes(x=elevation, y=delta.boom))  +
  xlab("") +
  ylab("") +
  my.theme

plot.call.elev <- ggplot(newdat.elev) +
  geom_ribbon(aes(x=elevation, ymin=delta.call.low, ymax = delta.call.high), alpha = 0.5) +
  geom_line(aes(x=elevation, y=delta.call)) +
  xlab("Elevation") +
  ylab("") +
  my.theme

ggsave(grid.arrange(plot.boom.grass, plot.boom.elev, plot.call.grass, plot.call.elev, ncol=2, nrow=2), 
       file = "Figs/CovariatePredictions.jpeg", width = 8, height = 8, device="jpeg")

#FIGURE 6. SPATIAL PREDICTIONS####
pred <- read.csv("SpatialPredictions.csv")

plot.density.boom <- ggplot(pred) +
  geom_raster(aes(x = x, y = y, fill=densityha.boom), na.rm=TRUE) +
  scale_fill_viridis_c(name="Territorial\nmales\nper ha") +
  xlab("") +
  ylab("") +
  my.theme +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot.density.call <- ggplot(pred) +
  geom_raster(aes(x = x, y = y, fill=densityha.call), na.rm=TRUE) +
  scale_fill_viridis_c(name="Total\nmales\nper ha") + 
  xlab("") +
  ylab("") +
  my.theme +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(grid.arrange(plot.density.boom, plot.density.call, ncol = 2), 
       file = "Figs/SpatialPredictions.jpeg", width = 12, height = 5, device="jpeg")


####### SUMMARY STATS FOR THINGS#####

#Number of surveys for each type
dat <- read.csv("SurveyDataWithCovs.csv") %>% 
  mutate(Kenow = ifelse(FireHistory==2017, "impact", "control"),
         Kenow = ifelse(is.na(Kenow), "control", Kenow)) %>% 
  dplyr::select(survey, station, year, X, Y) %>% 
  unique()
table(dat$survey)
table(dat$survey, dat$year)

#number of survey stations & years
dat.st <- dat %>% 
  dplyr::select(station, year, survey) %>% 
  unique() %>% 
  group_by(station, survey) %>% 
  summarize(n=n()) %>% 
  ungroup()

table(dat.st$survey)
table(dat.st$survey, dat.st$n)
  

#Number of detections
dat <- read.csv("SurveyDataWithCovs.csv")
  group_by(survey, station, year, X, Y) %>% 
  summarize(booms = sum(boom),
            calls = sum(call)) %>% 
  ungroup() %>% 
  mutate(boom.det = ifelse(booms > 0, 1, 0),
         call.det = ifelse(calls > 0, 1, 0))

table(dat$survey, dat$boom.det)    
table(dat$survey, dat$call.det)    
