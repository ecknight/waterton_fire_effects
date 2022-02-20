library(tidyverse)
library(raster)
library(spatialEco)
library(sf)
library(gridExtra)

#1. Read in model predictions----
newdat <- read.csv("OccupancyModelPredictions.csv")

#2. Read in spatial data----
grass.300 <- raster("rasters/grass-300.tif") %>% 
  aggregate(3)
names(grass.300) <- "grass"
elevation <- raster("rasters/DEM_10m.tif") %>% 
  aggregate(3)
names(elevation) <- "elevation"

#3. Predictions---
lambda <- read.csv("LambdaEstimates.csv")
edr <- read.csv("EDR.csv")[1,]$eda

pred <- grass.300 %>% 
  as.data.frame(xy=TRUE) %>% 
  cbind(elevation %>% 
          as.data.frame(xy=TRUE) %>% 
          dplyr::select(elevation)) %>% 
  dplyr::filter(!is.na(grass),
                !is.na(elevation),
                elevation > 0) %>% 
  mutate(grass = round(grass, 2),
         elevation = round(elevation)) %>% 
  left_join(newdat) %>% 
  mutate(lambda.boom = lambda$lambda[1],
         lambda.call = lambda$lambda[2],
         densityha.boom = delta.boom*lambda.boom/edr,
         densityha.call = delta.call*lambda.call/edr,
         cellsize = 900,
         densitycell.boom = densityha.boom*cellsize/10000,
         densitycell.call = densityha.call*cellsize/10000)
  
write.csv(pred, "SpatialPredictions.csv", row.names = FALSE)

#3. Population estimates----
pop <- pred %>% 
  summarize(pop.boom = sum(densitycell.boom, na.rm=TRUE),
            pop.call = sum(densitycell.call, na.rm=TRUE),
            cells= n())
pop

#4. Plot----
plot.density.boom <- ggplot(pred) +
  geom_raster(aes(x = x, y = y, fill=densityha.boom), na.rm=TRUE) +
  scale_fill_viridis_c(name="Territorial males\nper ha")
plot.density.boom

plot.density.call <- ggplot(pred) +
  geom_raster(aes(x = x, y = y, fill=densityha.call), na.rm=TRUE) +
  scale_fill_viridis_c(name="Total males\nper ha")
plot.density.call

ggsave(grid.arrange(plot.density.boom, plot.density.call, ncol = 2), 
       file = "Figs/SpatialPredictions.jpeg", width = 16, height = 7, device="jpeg")

pred.long <- pred %>% 
  pivot_longer(densityha.boom:densityha.call, names_to = "response", values_to = "density", names_prefix = "densityha.")

plot.density <- ggplot(pred.long) +
  geom_raster(aes(x=x, y=y, fill=density), na.rm=TRUE) +
  scale_fill_viridis_c(name="Males\nper ha") +
  facet_wrap(~response)

ggsave(plot.density, 
       file = "Figs/SpatialPredictions_Facet.jpeg", width = 16, height = 7, device="jpeg")

#5. Territory size----
boom.max <- max(pred$densitycell.boom)
territory <- 1/boom.max*0.4229*0.3041259
