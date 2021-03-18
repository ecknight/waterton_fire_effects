#title: Visualization of available WLNP CONI data
#author: Elly C. Knight
#created: March 17, 2021

library(tidyverse)
library(sf)
library(ggmap)
library(ggspatial)

aru <- read.csv("ARUData.csv") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  dplyr::select(SurveyType, Year, Site, Presence)

hum.11 <- read.csv("HumanData.csv") %>% 
  filter(Zone=="11U") %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26911) %>% 
  st_transform(crs=4326) %>% 
  dplyr::select(SurveyType, Year, Site, Presence)

hum.12 <- read.csv("HumanData.csv") %>% 
  filter(Zone=="12U") %>% 
  st_as_sf(coords=c("Easting", "Northing"), crs=26912) %>% 
  st_transform(crs=4326) %>% 
  dplyr::select(SurveyType, Year, Site, Presence) %>% 
  dplyr::select(SurveyType, Year, Site, Presence)

dat <- rbind(aru, hum.11, hum.12) %>% 
  st_coordinates() %>% 
  cbind(rbind(aru, hum.11, hum.12)) %>% 
  mutate(surveyID=paste0(Year,"-",SurveyType))

#Get background data
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

dat.map <- ggmap(map_transparent) +
  geom_spatial_point(aes(x = X, y = Y,
                         colour=factor(Presence)),
                     data = dat, 
                     crs=4326,
                     alpha = 0.7,
                     size=2,
                     show.legend = TRUE) +
  facet_wrap(~surveyID, nrow=1) +
  scale_colour_manual(values=c("blue", "orange"))

ggsave(dat.map, file="OccurrenceVisualization.jpeg", width=13, height=4)
