#title: ZIP abundance model fitting for CONI density modelling
#author: Elly C. Knight, adapted from Peter Solymos
#created: November 8, 2021

library(tidyverse)
library(lubridate)
library(AICcmodavg)

options(scipen=99999)

source("src/mvocc.R")

#1. Wrangle----
dat <- read.csv("SurveyDataWithOffsets.csv")  %>% 
  mutate(ID=paste0(survey,"-",station,"-",year),
         boom=ifelse(detection==2, 1, 0),
         call=ifelse(detection>0, 1, 0)) %>% 
  arrange(ID, DateTime) %>% 
  group_by(ID) %>% 
  mutate(n=row_number()) %>% 
  ungroup()

#2. Add covariates to data----
cov <- read.csv("SurveyLocationsWithCovs.csv") %>% 
  mutate(X=round(X),
         Y=round(Y)) %>% 
  unique()

dat.cov <- dat %>% 
  mutate(station = gsub(pattern="HORSE-0", replacement="HORSE-", x=station),
         station = gsub(pattern="LAKE-0", replacement="LAKE-", x=station),
         station = gsub(pattern="OIL-0", replacement="OIL-", x=station),
         station = gsub(pattern="RRP-0", replacement="RRP-", x=station),
         station = gsub(pattern="CLOUDY-0", replacement="CLOUDY-", x=station),
         station = gsub(pattern="WLNP-006-001", replacement="WLNP-6-1", x=station),
         station = gsub(pattern="WLNP-007-001", replacement="WLNP-7-1", x=station),
         station = gsub(pattern="CONI-00", replacement="CONI-", x=station),
         station = gsub(pattern="CONI-0", replacement="CONI-", x=station),
         station = gsub(pattern="-00", replacement="-", x=station),
         station = gsub(pattern="-0", replacement="-", x=station)) %>% 
  inner_join(cov) %>% 
  unique() %>% 
  mutate(FireTime=year-FireHistory,
         FireTime=ifelse(FireTime < 0, NA, FireTime),
         FireTime=ifelse(is.na(FireTime), 200, FireTime),
         FireSeverity = ifelse(year > 2017, FireSeverity, NA))  %>% 
  dplyr::filter(!is.na(pine.300),
                !is.na(sand.300),
                !is.na(cover.300))

write.csv(dat.cov, "SurveyDataWithOffsets&Covariates.csv", row.names = FALSE)

#3. Format for occupancy----
lambda <- read.csv("LambdaEstimates.csv")

#Boom
#dat.boom <- dat.cov %>% 
#  dplyr::filter(survey=="ARU")

dat.boom <- dat.cov

station.boom <- dat.boom %>% 
  dplyr::select(station, year, fire, ID, X, Y, Elevation, FireTime, grass.300, sand.300, trails.300, cover.300, Elevation) %>% 
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
dat.call <- dat.cov

station.call <- dat.call %>% 
  dplyr::select(station, year, fire, ID, X, Y, Elevation, FireTime, grass.300, sand.300, trails.300, cover.300, Elevation) %>% 
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

#4. Model boom----
#method <- "Nelder-Mead" # fail fast
method <- "SANN"
#method <- "DE" # slow and sure

X.boom <- matrix(1, M.boom, 1)
Xveg1.boom <- model.matrix(~grass.300*cover.300 + sand.300, station.boom)
Xveg2.boom <- model.matrix(~grass.300*cover.300, station.boom)
Xveg3.boom <- model.matrix(~grass.300 + cover.300 + sand.300, station.boom)
Xveg4.boom <- model.matrix(~grass.300 + cover.300, station.boom)
Xveg5.boom <- model.matrix(~grass.300 + sand.300, station.boom)
Xveg6.boom <- model.matrix(~cover.300 + sand.300, station.boom)
Xveg7.boom <- model.matrix(~grass.300, station.boom)
Xveg8.boom <- model.matrix(~cover.300, station.boom)
Xveg9.boom <- model.matrix(~sand.300, station.boom)

#Xfire.boom <- model.matrix(~FireTime, station.boom)

#Boom
mod.null.boom <- mvocc(y.boom, X.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg1.boom <- mvocc(y.boom, Xveg1.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg2.boom <- mvocc(y.boom, Xveg2.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg3.boom <- mvocc(y.boom, Xveg3.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg4.boom <- mvocc(y.boom, Xveg4.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg5.boom <- mvocc(y.boom, Xveg5.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg6.boom <- mvocc(y.boom, Xveg6.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg7.boom <- mvocc(y.boom, Xveg7.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg8.boom <- mvocc(y.boom, Xveg8.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg9.boom <- mvocc(y.boom, Xveg9.boom, Z.boom, p.boom, lamvec.boom, method=method)

aic <- AIC(mod.null.boom, mod.veg1.boom, mod.veg2.boom, mod.veg3.boom, mod.veg4.boom, mod.veg5.boom, mod.veg6.boom, mod.veg7.boom, mod.veg8.boom, mod.veg9.boom)
aic

mod.veg.boom <- mod.veg5.boom

Xvegelev.boom <- model.matrix(~grass.300 + sand.300 + Elevation, station.boom)
Xvegtrails.boom <- model.matrix(~grass.300 + sand.300 + trails.300, station.boom)

mod.vegelev.boom <- mvocc(y.boom, Xvegelev.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.vegtrails.boom <- mvocc(y.boom, Xvegtrails.boom, Z.boom, p.boom, lamvec.boom, method=method)

aic <- AIC(mod.veg.boom, mod.vegelev.boom, mod.vegtrails.boom)
aic

Xvegfire.boom <- model.matrix(~grass.300 + sand.300 + FireTime, station.boom)

mod.vegfire.boom <- mvocc(y.boom, Xvegfire.boom, Z.boom, p.boom, lamvec.boom, method=method)

aic <- AIC(mod.veg.boom, mod.vegfire.boom)
aic

best.boom <- mod.veg.boom
summary(best.boom)

saveRDS(best.boom, "OccupancyModel_Boom.rds")

#5. Model call----
#method <- "Nelder-Mead" # fail fast
method <- "SANN"
#method <- "DE" # slow and sure

X.call <- matrix(1, M.call, 1)
Xveg1.call <- model.matrix(~grass.300*cover.300 + sand.300, station.call)
Xveg2.call <- model.matrix(~grass.300*cover.300, station.call)
Xveg3.call <- model.matrix(~grass.300 + cover.300 + sand.300, station.call)
Xveg4.call <- model.matrix(~grass.300 + cover.300, station.call)
Xveg5.call <- model.matrix(~grass.300 + sand.300, station.call)
Xveg6.call <- model.matrix(~cover.300 + sand.300, station.call)
Xveg7.call <- model.matrix(~grass.300, station.call)
Xveg8.call <- model.matrix(~cover.300, station.call)
Xveg9.call <- model.matrix(~sand.300, station.call)

#Xfire.call <- model.matrix(~FireTime, station.call)

#call
mod.null.call <- mvocc(y.call, X.call, Z.call, p.call, lamvec.call, method=method)
mod.veg1.call <- mvocc(y.call, Xveg1.call, Z.call, p.call, lamvec.call, method=method)
mod.veg2.call <- mvocc(y.call, Xveg2.call, Z.call, p.call, lamvec.call, method=method)
mod.veg3.call <- mvocc(y.call, Xveg3.call, Z.call, p.call, lamvec.call, method=method)
mod.veg4.call <- mvocc(y.call, Xveg4.call, Z.call, p.call, lamvec.call, method=method)
mod.veg5.call <- mvocc(y.call, Xveg5.call, Z.call, p.call, lamvec.call, method=method)
mod.veg6.call <- mvocc(y.call, Xveg6.call, Z.call, p.call, lamvec.call, method=method)
mod.veg7.call <- mvocc(y.call, Xveg7.call, Z.call, p.call, lamvec.call, method=method)
mod.veg8.call <- mvocc(y.call, Xveg8.call, Z.call, p.call, lamvec.call, method=method)
mod.veg9.call <- mvocc(y.call, Xveg9.call, Z.call, p.call, lamvec.call, method=method)

aic <- AIC(mod.null.call, mod.veg1.call, mod.veg2.call, mod.veg3.call, mod.veg4.call, mod.veg5.call, mod.veg6.call, mod.veg7.call, mod.veg8.call, mod.veg9.call)
aic

mod.veg.call <- mod.veg7.call

Xvegelev.call <- model.matrix(~grass.300 + Elevation, station.call)
Xvegtrails.call <- model.matrix(~grass.300 + trails.300, station.call)

mod.vegelev.call <- mvocc(y.call, Xvegelev.call, Z.call, p.call, lamvec.call, method=method)
mod.vegtrails.call <- mvocc(y.call, Xvegtrails.call, Z.call, p.call, lamvec.call, method=method)

aic <- AIC(mod.veg.call, mod.vegelev.call, mod.vegtrails.call)
aic

Xvegfire.call <- model.matrix(~grass.300 + FireTime, station.call)

mod.vegfire.call <- mvocc(y.call, Xvegfire.call, Z.call, p.call, lamvec.call, method=method)

aic <- AIC(mod.veg.call, mod.vegfire.call)
aic

best.call <- mod.veg.call
summary(best.call)

saveRDS(best.call, "OccupancyModel_Call.rds")

#6. Estimates----
best.boom$coef
best.call$coef

#6a. Suitability----
delta.boom <- data.frame(delta = plogis(drop(Xveg5.boom %*% best.boom$coef[1:3])),
                         sa=station.boom$station)
summary(delta.boom$delta)

delta.call <- data.frame(delta = plogis(drop(Xveg7.call %*% best.call$coef[1:2])),
                         sa=station.call$station)
summary(delta.call$delta)

#6b. Occupied----
summary(1-exp(-lambda$lambda[1]))
summary(1-exp(-lambda$lambda[2]))

#6c. Active----
summary(as.numeric(p.boom))
summary(as.numeric(p.call))

#6d. Used----
plogis(best.boom$coef[4])
plogis(best.call$coef[3])

#6e. Density----
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
