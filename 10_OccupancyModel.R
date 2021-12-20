#title: Modified occupancy model fitting for CONI density modelling
#author: Elly C. Knight, adapted from Peter Solymos
#created: November 8, 2021

library(tidyverse)
library(lubridate)
#library(tidylog)
library(AICcmodavg)
library(sf)
library(gridExtra)

options(scipen=99999)

source("src/mvocc.R")

#1. Wrangle----
dat.raw <- read.csv("SurveyDataWithCovs.csv") %>% 
  arrange(ID, DateTime) %>% 
  group_by(ID) %>% 
  mutate(n=row_number()) %>% 
  ungroup()

dat.occ <- dat.raw %>% 
  group_by(ID, survey) %>% 
  summarize(occ.boom = ifelse(sum(boom) > 0, 1, 0),
            occ.call = ifelse(sum(call) > 0, 1, 0)) %>% 
  ungroup()

table(dat.occ$occ.boom, dat.occ$occ.call, dat.occ$survey)
table(dat.occ$survey)

dat <- dat.raw %>% 
  left_join(dat.occ) %>% 
  mutate(elevation.s = scale(Elevation),
         grass.s = scale(grass.300),
         cover.s = scale(cover.300),
         sand.s = scale(sand.300),
         wetland.s = scale(wetland.300),
         firetime.s = scale(FireTime))

#4. Visualize----
ggplot(dat, aes(x=grass.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=grass.300, y=call)) +
  geom_smooth()

ggplot(dat, aes(x=cover.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=cover.300, y=call)) +
  geom_smooth()

ggplot(dat, aes(x=sand.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=sand.300, y=call)) +
  geom_smooth()

ggplot(dat, aes(x=wetland.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=wetland.300, y=call)) +
  geom_smooth()

ggplot(dat, aes(x=Elevation, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=Elevation, y=call)) +
  geom_smooth()

ggplot(dat, aes(x=FireTime, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=FireTime, y=call)) +
  geom_smooth()


#4. Format for occupancy----
lambda <- read.csv("LambdaEstimates.csv")

#Boom
dat.boom <- dat %>% 
  dplyr::filter(survey=="ARU")

station.boom <- dat.boom %>% 
  dplyr::select(station, year, fire, ID, X, Y, firetime.s, grass.s, elevation.s) %>% 
  unique()

y.boom <- dat.boom %>% 
  dplyr::select(ID, n, boom) %>% 
  spread(key = n, value = boom) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.matrix()

p.boom <- dat.boom %>% 
  dplyr::select(ID, n, p.boom) %>% 
  spread(key = n, value = p.boom) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.matrix()

M.boom <- nrow(y.boom)
J.boom <- ncol(y.boom)
Z.boom <- matrix(1, M.boom*J.boom, 1)
lamvec.boom <- rep(lambda$lambda[1], nrow(station.boom))

#Call
dat.call <- dat

station.call <- dat.call %>% 
  dplyr::select(station, year, fire, ID, X, Y, firetime.s, grass.s, elevation.s) %>% 
  unique()

y.call <- dat.call %>% 
  dplyr::select(ID, n, call) %>% 
  spread(key = n, value = call) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.matrix()

p.call <- dat.call %>% 
  dplyr::select(ID, n, p.call) %>% 
  spread(key = n, value = p.call) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.matrix()

M.call <- nrow(y.call)
J.call <- ncol(y.call)
Z.call <- matrix(1, M.call*J.call, 1)
lamvec.call <- rep(lambda$lambda[1], nrow(station.call))

#6. Confirm polynomial for elevation----
#method <- "Nelder-Mead" # fail fast
#method <- "SANN"
method <- "DE" # slow and sure

Xelev1.boom <- model.matrix(~elevation.s, station.boom)
Xelev2.boom <- model.matrix(~poly(elevation.s, 2), station.boom)

mod.elev1.boom <- mvocc(y.boom, Xelev1.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.elev2.boom <- mvocc(y.boom, Xelev2.boom, Z.boom, p.boom, lamvec.boom, method=method)

AIC(mod.elev1.boom, mod.elev2.boom)

Xelev1.call <- model.matrix(~elevation.s, station.call)
Xelev2.call <- model.matrix(~poly(elevation.s, 2), station.call)

mod.elev1.call <- mvocc(y.call, Xelev1.call, Z.call, p.call, lamvec.call, method=method)
mod.elev2.call <- mvocc(y.call, Xelev2.call, Z.call, p.call, lamvec.call, method=method)

AIC(mod.elev1.call, mod.elev2.call)

#7. Model matrices----
X.boom <- matrix(1, M.boom, 1)
Xveg1.boom <- model.matrix(~grass.s + elevation.s + firetime.s, station.boom)
Xveg2.boom <- model.matrix(~grass.s + firetime.s, station.boom)
Xveg3.boom <- model.matrix(~elevation.s + firetime.s, station.boom)
Xveg4.boom <- model.matrix(~grass.s + elevation.s, station.boom)
Xveg5.boom <- model.matrix(~grass.s, station.boom)
Xveg6.boom <- model.matrix(~firetime.s, station.boom)
Xveg7.boom <- model.matrix(~elevation.s, station.boom)

X.call <- matrix(1, M.call, 1)
Xveg1.call <- model.matrix(~grass.s + elevation.s + firetime.s, station.call)
Xveg2.call <- model.matrix(~grass.s + firetime.s, station.call)
Xveg3.call <- model.matrix(~elevation.s + firetime.s, station.call)
Xveg4.call <- model.matrix(~grass.s + elevation.s, station.call)
Xveg5.call <- model.matrix(~grass.s, station.call)
Xveg6.call <- model.matrix(~firetime.s, station.call)
Xveg7.call <- model.matrix(~elevation.s, station.call)

#7. Run models----
#method <- "Nelder-Mead" # fail fast
#method <- "SANN"
method <- "DE" # slow and sure

#Boom
set.seed(1234)
mod.null.boom <- mvocc(y.boom, X.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg1.boom <- mvocc(y.boom, Xveg1.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg2.boom <- mvocc(y.boom, Xveg2.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg3.boom <- mvocc(y.boom, Xveg3.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg4.boom <- mvocc(y.boom, Xveg4.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg5.boom <- mvocc(y.boom, Xveg5.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg6.boom <- mvocc(y.boom, Xveg6.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg7.boom <- mvocc(y.boom, Xveg7.boom, Z.boom, p.boom, lamvec.boom, method=method)

set.seed(1234)
mod.null.call <- mvocc(y.call, X.call, Z.call, p.call, lamvec.call, method=method)
mod.veg1.call <- mvocc(y.call, Xveg1.call, Z.call, p.call, lamvec.call, method=method)
mod.veg2.call <- mvocc(y.call, Xveg2.call, Z.call, p.call, lamvec.call, method=method)
mod.veg3.call <- mvocc(y.call, Xveg3.call, Z.call, p.call, lamvec.call, method=method)
mod.veg4.call <- mvocc(y.call, Xveg4.call, Z.call, p.call, lamvec.call, method=method)
mod.veg5.call <- mvocc(y.call, Xveg5.call, Z.call, p.call, lamvec.call, method=method)
mod.veg6.call <- mvocc(y.call, Xveg6.call, Z.call, p.call, lamvec.call, method=method)
mod.veg7.call <- mvocc(y.call, Xveg7.call, Z.call, p.call, lamvec.call, method=method)

#8. Pick best model----
mods.boom <- list(mod.null.boom, mod.veg1.boom, mod.veg2.boom, mod.veg3.boom, mod.veg4.boom, mod.veg5.boom, mod.veg6.boom, mod.veg7.boom)

aic.boom <- data.frame(df=sapply(mods.boom, function(z) length(coef(z))),
                       AIC=sapply(mods.boom, AIC),
                       loglik = sapply(mods.boom, function(z) logLik(z))) %>% 
  mutate(AICc = AIC + (2*df^2+2*df) / (nobs(mod.null.boom)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight) %>% 
  mutate(response = "boom")
aic.boom$model <- row.names(aic.boom)
aic.boom$model <- c("mod.null.boom", "mod.veg1.boom", "mod.veg2.boom", "mod.veg3.boom", "mod.veg4.boom", "mod.veg5.boom", "mod.veg6.boom", "mod.veg7.boom")

mods.call <- list(mod.null.call, mod.veg1.call, mod.veg2.call, mod.veg3.call, mod.veg4.call, mod.veg5.call, mod.veg6.call, mod.veg7.call)

aic.call <- data.frame(df=sapply(mods.call, function(z) length(coef(z))),
                       AIC=sapply(mods.call, AIC),
                       loglik = sapply(mods.call, function(z) logLik(z))) %>% 
  mutate(AICc = AIC + (2*df^2+2*df) / (nobs(mod.null.call)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight) %>% 
  mutate(response = "call")
aic.call$model <- c("mod.null.call", "mod.veg1.call", "mod.veg2.call", "mod.veg3.call", "mod.veg4.call", "mod.veg5.call", "mod.veg6.call", "mod.veg7.call")

aic <- rbind(aic.boom, aic.call) %>% 
  arrange(response, delta)
View(aic)

write.csv(aic, "OccupancyModelSelection.csv", row.names = FALSE)

summary(mod.veg5.boom)
best.boom <- mod.veg5.boom
saveRDS(best.boom, "OccupancyModel_Boom.rds")

summary(mod.veg4.call)
best.call <- mod.veg4.call
saveRDS(best.call, "OccupancyModel_Call.rds")

save.image("DEOccupancyModels.Rdata")
load("DEOccupancyModels.Rdata")

#9. Model predictions----
#Newdata
newdat <- data.frame(expand_grid(grass = seq(0, 1, 0.01),
                                 elevation = seq(12, 2938, 1))) %>% 
  mutate(grass.s = (grass - mean(dat$grass.300))/sd(dat$grass.300),
         elevation.s = (elevation - mean(dat$Elevation))/sd(dat$Elevation))

#Predict
Xnewdat.boom <- model.matrix(~grass.s, newdat)
newdat$delta.boom <- plogis(drop(Xnewdat.boom %*% best.boom$coef[1:2]))
#newdat.boom$delta.se <- plogis(drop(Xnewdat.boom %*% best.boom[["summary"]][7:8]))

Xnewdat.call <- model.matrix(~grass.s + elevation.s, newdat)
newdat$delta.call <- plogis(drop(Xnewdat.call %*% best.call$coef[1:3]))

write.csv(newdat, "OccupancyModelPredictions.csv", row.names = FALSE)

#Plot
newdat.grass <- newdat %>% 
  dplyr::filter(elevation==1500)

newdat.elev <- newdat %>% 
  dplyr::filter(grass==0.5)

plot.boom.grass <- ggplot(newdat.grass) +
  geom_line(aes(x=grass, y=delta.boom))

plot.call.grass <- ggplot(newdat.grass) +
  geom_line(aes(x=grass, y=delta.call))

plot.call.elev <- ggplot(newdat.elev) +
  geom_line(aes(x=elevation, y=delta.call))

grid.arrange(plot.boom.grass, plot.call.grass, plot.call.elev)


#10. Estimates----
best.boom$coef
best.call$coef

#10a. Suitability----
delta.boom <- data.frame(delta = plogis(drop(Xveg5.boom %*% best.boom$coef[1:2])),
                         sa=station.boom$station)
summary(delta.boom$delta)

delta.call <- data.frame(delta = plogis(drop(Xveg4.call %*% best.call$coef[1:3])),
                         sa=station.call$station)
summary(delta.call$delta)

#10b. Occupied----
summary(1-exp(-lambda$lambda[1]))
summary(1-exp(-lambda$lambda[2]))

#10c. Active----
summary(as.numeric(p.boom))
summary(as.numeric(p.call))

#10d. Used----
plogis(best.boom$coef[3])
plogis(best.call$coef[4])

#10e. Density----
edr <- read.csv("EDR.csv")[1,]$eda

estimates.boom <- delta.boom %>% 
  cbind(lamvec.boom) %>% 
  rename(lambda = lamvec.boom) %>% 
  mutate(density = delta*lambda/edr)
mean(estimates.boom$density)

estimates.call <- delta.call %>% 
  cbind(lamvec.call) %>% 
  rename(lambda = lamvec.call) %>% 
  mutate(density = delta*lambda/edr)
mean(estimates.call$density)
