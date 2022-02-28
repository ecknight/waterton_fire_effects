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
         lambda.boom.high = lambda$quanthigh[1],
         lambda.boom.low = lambda$quantlow[1],
         lambda.call = lambda$lambda[2],
         lambda.call.high = lambda$quanthigh[2],
         lambda.call.low = lambda$quantlow[2],
         densityha.boom = delta.boom*lambda.boom/edr,
         densityha.boom.high = delta.boom.high*lambda.boom.high/edr,
         densityha.boom.low = delta.boom.low*lambda.boom.low/edr,
         densityha.call = delta.call*lambda.call/edr,
         densityha.call.high = delta.call.high*lambda.call.high/edr,
         densityha.call.low = delta.call.low*lambda.call.low/edr,
         cellsize = 900,
         densitycell.boom = densityha.boom*cellsize/10000,
         densitycell.boom.high = densityha.boom.high*cellsize/10000,
         densitycell.boom.low = densityha.boom.low*cellsize/10000,
         densitycell.call = densityha.call*cellsize/10000,
         densitycell.call.high = densityha.call.high*cellsize/10000,
         densitycell.call.low = densityha.call.low*cellsize/10000)
  
write.csv(pred, "SpatialPredictions.csv", row.names = FALSE)

#3. Population estimates----
pop <- pred %>% 
  summarize(pop.boom = sum(densitycell.boom, na.rm=TRUE),
            pop.boom.high = sum(densitycell.boom.high, na.rm=TRUE),
            pop.boom.low = sum(densitycell.boom.low, na.rm=TRUE),
            pop.call = sum(densitycell.call, na.rm=TRUE),
            pop.call.high = sum(densitycell.call.high, na.rm=TRUE),
            pop.call.low = sum(densitycell.call.low, na.rm=TRUE),
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
#For most dense areas
territory <- 1/max(pred$densityha.boom)
territory.high <- 1/max(pred$densityha.boom.high)
territory.low <- 1/max(pred$densityha.boom.low)
