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
         FireSeverity = ifelse(year > 2017, FireSeverity, NA)) 

#3. Format for occupancy----
lambda <- read.csv("LambdaEstimates.csv")

#Boom first
dat.boom <- dat.cov %>% 
  dplyr::filter(survey=="ARU")

station.boom <- dat.boom %>% 
  dplyr::select(station, year, fire, ID, X, Y, Elevation, FireHistory, FireSeverity, Soil, Vegetation, VegetationCover, VegetationHeight, FireTime) %>% 
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

X.boom <- matrix(1, M.boom, 1)
Xfire.boom <- model.matrix(~FireTime, station.boom)

#Call
dat.call <- dat.cov

station.call <- dat.call %>% 
  dplyr::select(station, year, fire, ID, X, Y, Elevation, FireHistory, FireSeverity, Soil, Vegetation, VegetationCover, VegetationHeight, FireTime) %>% 
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
lamvec.call <- rep(lambda$lambda[2], nrow(station.call))

X.call <- matrix(1, M.call, 1)
Xfire.call <- model.matrix(~FireTime, station.call)

#4. Model----
#method <- "Nelder-Mead" # fail fast
#method <- "SANN"
method <- "DE" # slow and sure

#Boom
o00.boom <- mvocc(y.boom, X.boom, Z.boom, p.boom, lamvec.boom, method=method)
ofire.boom <- mvocc(y.boom, Xfire.boom, Z.boom, p.boom, lamvec.boom, method=method)

aic <- AIC(o00.boom, ofire.boom)
aic

best.boom <- o00.boom

#Call
o00.call <- mvocc(y.call, X.call, Z.call, p.call, lamvec.call, method=method)
ofire.call <- mvocc(y.call, Xfire.call, Z.call, p.call, lamvec.call, method=method)

aic <- AIC(o00.call, ofire.call)
aic

best.call <- o00.call

#5. Probability of used----
plogis(best.boom$coef[max(length(best.boom$coef))])
plogis(best.call$coef[max(length(best.call$coef))]) # used
