library(tidyverse)
library(raster)
library(spatialEco)
library(sf)
library(gridExtra)

best.boom <- readRDS("OccupancyModel_Boom.rds")
best.call <- readRDS("OccupancyModel_Call.rds")

#1. Predictions for new variable combos----
newdat.boom <- data.frame(expand_grid(grass.300=seq(0, 1, 0.01),
                                 sand.300 = seq(0, 1, 0.01)))
Xnewdat.boom <- model.matrix(~grass.300 + sand.300, newdat.boom)
newdat.boom$delta <- plogis(drop(Xnewdat.boom %*% best.boom$coef[1:3]))

newdat.call <- data.frame(expand_grid(grass.300=seq(0, 1, 0.01)))
Xnewdat.call <- model.matrix(~grass.300, newdat.call)
newdat.call$delta <- plogis(drop(Xnewdat.call %*% best.call$coef[1:2]))

#2. Read in spatial data----
sand.300 <- raster("rasters/sand-300.tif")
names(sand.300) <- "sand.300"
grass.300 <- raster("rasters/grass-300.tif")
names(grass.300) <- "grass.300"

#3. Predictions---
lambda <- read.csv("LambdaEstimates.csv")
edr <- read.csv("EDR.csv")[1,]$eda

pred <- sand.300 %>% 
  as.data.frame(xy=TRUE) %>% 
  cbind(grass.300 %>% 
          as.data.frame(xy=TRUE) %>% 
          dplyr::select(grass.300)) %>% 
  dplyr::filter(!is.na(grass.300)) %>% 
  mutate(grass.300 = round(grass.300, 2),
         sand.300 = round(sand.300, 2)) %>% 
  left_join(newdat.boom) %>% 
  rename(delta.boom = delta) %>% 
  left_join(newdat.call) %>% 
  rename(delta.call = delta) %>% 
  mutate(lambda.boom = lambda$lambda[1],
         lambda.call = lambda$lambda[2],
         densityha.boom = delta.boom*lambda.boom/edr,
         densityha.call = delta.call*lambda.call/edr,
         cellsize = 100,
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
