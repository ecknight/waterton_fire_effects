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
  mutate(n=row_number(),
         extra = detection==3) %>% 
  ungroup()

dat.occ <- dat.raw %>% 
  group_by(ID, survey) %>% 
  summarize(occ.boom = ifelse(sum(boom) > 0, 1, 0),
            occ.call = ifelse(sum(call) > 0, 1, 0),
            occ.extra = ifelse(occ.boom==0 & occ.call==1, 1, 0)) %>% 
  ungroup()

table(dat.occ$occ.boom, dat.occ$occ.call, dat.occ$survey)
table(dat.occ$occ.boom, dat.occ$occ.extra, dat.occ$survey)
table(dat.occ$survey)

dat <- dat.raw %>% 
  left_join(dat.occ) %>% 
  mutate(elevation.s = scale(Elevation),
         grass.s = scale(grass.300),
         cover.s = scale(cover.300),
         sand.s = scale(sand.300),
         wetland.s = scale(wetland.300),
         firetime.s = scale(FireTime),
         evi.s = scale(evi)) %>% 
  mutate(extra = ifelse(extra==FALSE, 0, 1))

#4. Visualize----
ggplot(dat, aes(x=grass.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=grass.300, y=call)) +
  geom_smooth()
ggplot(dat, aes(x=grass.300, y=extra)) +
  geom_smooth()

ggplot(dat, aes(x=cover.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=cover.300, y=call)) +
  geom_smooth()
ggplot(dat, aes(x=cover.300, y=extra)) +
  geom_smooth()

ggplot(dat, aes(x=sand.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=sand.300, y=call)) +
  geom_smooth()
ggplot(dat, aes(x=sand.300, y=extra)) +
  geom_smooth()

ggplot(dat, aes(x=wetland.300, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=wetland.300, y=call)) +
  geom_smooth()
ggplot(dat, aes(x=wetland.300, y=extra)) +
  geom_smooth()

ggplot(dat, aes(x=Elevation, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=Elevation, y=call)) +
  geom_smooth()
ggplot(dat, aes(x=Elevation, y=extra)) +
  geom_smooth()

ggplot(dat, aes(x=FireTime, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=FireTime, y=call)) +
  geom_smooth()
ggplot(dat, aes(x=FireTime, y=extra)) +
  geom_smooth()

ggplot(dat, aes(x=evi, y=boom)) +
  geom_smooth()
ggplot(dat, aes(x=evi, y=call)) +
  geom_smooth()
ggplot(dat, aes(x=evi, y=extra)) +
  geom_smooth()


#4. Format for occupancy----
lambda <- read.csv("LambdaEstimates.csv")

#Boom
dat.boom <- dat %>% 
  dplyr::filter(survey=="ARU")

station.boom <- dat.boom %>% 
  dplyr::select(station, year, fire, ID, X, Y, firetime.s, grass.s, elevation.s, wetland.s, evi.s, cover.s) %>% 
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
  dplyr::select(station, year, fire, ID, X, Y, firetime.s, grass.s, elevation.s, wetland.s, evi.s, cover.s) %>% 
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

#Extraterritorial
dat.extra <- dat %>% 
  dplyr::filter(survey=="ARU")

station.extra <- dat.extra %>% 
  dplyr::select(station, year, fire, ID, X, Y, firetime.s, grass.s, elevation.s, wetland.s, evi.s, cover.s) %>% 
  unique()

y.extra <- dat.extra %>% 
  dplyr::select(ID, n, extra) %>% 
  spread(key = n, value = extra) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.matrix()

p.extra <- dat.extra %>% 
  dplyr::select(ID, n, p.extra) %>% 
  spread(key = n, value = p.extra) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.matrix()

M.extra <- nrow(y.extra)
J.extra <- ncol(y.extra)
Z.extra <- matrix(1, M.extra*J.extra, 1)
lamvec.extra <- rep(lambda$lambda[3], nrow(station.extra))

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

Xelev1.extra <- model.matrix(~elevation.s, station.extra)
Xelev2.extra <- model.matrix(~poly(elevation.s, 2), station.extra)

mod.elev1.extra <- mvocc(y.extra, Xelev1.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.elev2.extra <- mvocc(y.extra, Xelev2.extra, Z.extra, p.extra, lamvec.extra, method=method)

AIC(mod.elev1.extra, mod.elev2.extra)

#7. Model matrices----
X.boom <- matrix(1, M.boom, 1)
Xveg1.boom <- model.matrix(~grass.s*evi.s + elevation.s + firetime.s, station.boom)
Xveg2.boom <- model.matrix(~grass.s + evi.s + elevation.s + firetime.s, station.boom)
Xveg3.boom <- model.matrix(~grass.s*evi.s + firetime.s, station.boom)
Xveg4.boom <- model.matrix(~grass.s*evi.s + elevation.s, station.boom)
Xveg5.boom <- model.matrix(~grass.s + evi.s + firetime.s, station.boom)
Xveg6.boom <- model.matrix(~grass.s + evi.s + elevation.s, station.boom)
Xveg7.boom <- model.matrix(~grass.s*evi.s, station.boom)
Xveg8.boom <- model.matrix(~grass.s + evi.s, station.boom)
Xveg9.boom <- model.matrix(~grass.s + elevation.s + firetime.s, station.boom)
Xveg10.boom <- model.matrix(~evi.s + elevation.s + firetime.s, station.boom)
Xveg11.boom <- model.matrix(~grass.s + elevation.s, station.boom)
Xveg12.boom <- model.matrix(~grass.s + firetime.s, station.boom)
Xveg13.boom <- model.matrix(~elevation.s + firetime.s, station.boom)
Xveg14.boom <- model.matrix(~evi.s + firetime.s, station.boom)
Xveg15.boom <- model.matrix(~elevation.s + evi.s, station.boom)
Xveg16.boom <- model.matrix(~grass.s, station.boom)
Xveg17.boom <- model.matrix(~evi.s, station.boom)
Xveg18.boom <- model.matrix(~elevation.s, station.boom)
Xveg19.boom <- model.matrix(~firetime.s, station.boom)

X.call <- matrix(1, M.call, 1)
Xveg1.call <- model.matrix(~grass.s*evi.s + elevation.s + firetime.s, station.call)
Xveg2.call <- model.matrix(~grass.s + evi.s + elevation.s + firetime.s, station.call)
Xveg3.call <- model.matrix(~grass.s*evi.s + firetime.s, station.call)
Xveg4.call <- model.matrix(~grass.s*evi.s + elevation.s, station.call)
Xveg5.call <- model.matrix(~grass.s + evi.s + firetime.s, station.call)
Xveg6.call <- model.matrix(~grass.s + evi.s + elevation.s, station.call)
Xveg7.call <- model.matrix(~grass.s*evi.s, station.call)
Xveg8.call <- model.matrix(~grass.s + evi.s, station.call)
Xveg9.call <- model.matrix(~grass.s + elevation.s + firetime.s, station.call)
Xveg10.call <- model.matrix(~evi.s + elevation.s + firetime.s, station.call)
Xveg11.call <- model.matrix(~grass.s + elevation.s, station.call)
Xveg12.call <- model.matrix(~grass.s + firetime.s, station.call)
Xveg13.call <- model.matrix(~elevation.s + firetime.s, station.call)
Xveg14.call <- model.matrix(~evi.s + firetime.s, station.call)
Xveg15.call <- model.matrix(~elevation.s + evi.s, station.call)
Xveg16.call <- model.matrix(~grass.s, station.call)
Xveg17.call <- model.matrix(~evi.s, station.call)
Xveg18.call <- model.matrix(~elevation.s, station.call)
Xveg19.call <- model.matrix(~firetime.s, station.call)

X.extra <- matrix(1, M.extra, 1)
Xveg1.extra <- model.matrix(~grass.s*cover.s + elevation.s + firetime.s, station.extra)
Xveg2.extra <- model.matrix(~grass.s + cover.s + elevation.s + firetime.s, station.extra)
Xveg3.extra <- model.matrix(~grass.s*cover.s + firetime.s, station.extra)
Xveg4.extra <- model.matrix(~grass.s*cover.s + elevation.s, station.extra)
Xveg5.extra <- model.matrix(~grass.s + cover.s + firetime.s, station.extra)
Xveg6.extra <- model.matrix(~grass.s + cover.s + elevation.s, station.extra)
Xveg7.extra <- model.matrix(~grass.s*cover.s, station.extra)
Xveg8.extra <- model.matrix(~grass.s + cover.s, station.extra)
Xveg9.extra <- model.matrix(~grass.s + elevation.s + firetime.s, station.extra)
Xveg10.extra <- model.matrix(~cover.s + elevation.s + firetime.s, station.extra)
Xveg11.extra <- model.matrix(~grass.s + elevation.s, station.extra)
Xveg12.extra <- model.matrix(~grass.s + firetime.s, station.extra)
Xveg13.extra <- model.matrix(~elevation.s + firetime.s, station.extra)
Xveg14.extra <- model.matrix(~cover.s + firetime.s, station.extra)
Xveg15.extra <- model.matrix(~elevation.s + cover.s, station.extra)
Xveg16.extra <- model.matrix(~grass.s, station.extra)
Xveg17.extra <- model.matrix(~cover.s, station.extra)
Xveg18.extra <- model.matrix(~elevation.s, station.extra)
Xveg19.extra <- model.matrix(~firetime.s, station.extra)

#7. Run models----
#method <- "Nelder-Mead" # fail fast
#method <- "SANN"
method <- "DE" # slow and sure

#Boom
set.seed(999)
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
mod.veg10.boom <- mvocc(y.boom, Xveg10.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg11.boom <- mvocc(y.boom, Xveg11.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg12.boom <- mvocc(y.boom, Xveg12.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg13.boom <- mvocc(y.boom, Xveg13.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg14.boom <- mvocc(y.boom, Xveg14.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg15.boom <- mvocc(y.boom, Xveg15.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg16.boom <- mvocc(y.boom, Xveg16.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg17.boom <- mvocc(y.boom, Xveg17.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg18.boom <- mvocc(y.boom, Xveg18.boom, Z.boom, p.boom, lamvec.boom, method=method)
mod.veg19.boom <- mvocc(y.boom, Xveg19.boom, Z.boom, p.boom, lamvec.boom, method=method)

set.seed(999)
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
mod.veg10.call <- mvocc(y.call, Xveg10.call, Z.call, p.call, lamvec.call, method=method)
mod.veg11.call <- mvocc(y.call, Xveg11.call, Z.call, p.call, lamvec.call, method=method)
mod.veg12.call <- mvocc(y.call, Xveg12.call, Z.call, p.call, lamvec.call, method=method)
mod.veg13.call <- mvocc(y.call, Xveg13.call, Z.call, p.call, lamvec.call, method=method)
mod.veg14.call <- mvocc(y.call, Xveg14.call, Z.call, p.call, lamvec.call, method=method)
mod.veg15.call <- mvocc(y.call, Xveg15.call, Z.call, p.call, lamvec.call, method=method)
mod.veg16.call <- mvocc(y.call, Xveg16.call, Z.call, p.call, lamvec.call, method=method)
mod.veg17.call <- mvocc(y.call, Xveg17.call, Z.call, p.call, lamvec.call, method=method)
mod.veg18.call <- mvocc(y.call, Xveg18.call, Z.call, p.call, lamvec.call, method=method)
mod.veg19.call <- mvocc(y.call, Xveg19.call, Z.call, p.call, lamvec.call, method=method)

set.seed(999)
mod.null.extra <- mvocc(y.extra, X.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg1.extra <- mvocc(y.extra, Xveg1.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg2.extra <- mvocc(y.extra, Xveg2.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg3.extra <- mvocc(y.extra, Xveg3.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg4.extra <- mvocc(y.extra, Xveg4.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg5.extra <- mvocc(y.extra, Xveg5.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg6.extra <- mvocc(y.extra, Xveg6.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg7.extra <- mvocc(y.extra, Xveg7.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg8.extra <- mvocc(y.extra, Xveg8.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg9.extra <- mvocc(y.extra, Xveg9.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg10.extra <- mvocc(y.extra, Xveg10.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg11.extra <- mvocc(y.extra, Xveg11.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg12.extra <- mvocc(y.extra, Xveg12.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg13.extra <- mvocc(y.extra, Xveg13.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg14.extra <- mvocc(y.extra, Xveg14.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg15.extra <- mvocc(y.extra, Xveg15.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg16.extra <- mvocc(y.extra, Xveg16.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg17.extra <- mvocc(y.extra, Xveg17.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg18.extra <- mvocc(y.extra, Xveg18.extra, Z.extra, p.extra, lamvec.extra, method=method)
mod.veg19.extra <- mvocc(y.extra, Xveg19.extra, Z.extra, p.extra, lamvec.extra, method=method)


#8. Pick best model----
mods.boom <- list(mod.null.boom, mod.veg1.boom, mod.veg2.boom, mod.veg3.boom, mod.veg4.boom, mod.veg5.boom, mod.veg6.boom, mod.veg7.boom, mod.veg8.boom, mod.veg9.boom, mod.veg10.boom, mod.veg11.boom, mod.veg12.boom, mod.veg13.boom, mod.veg14.boom, mod.veg15.boom, mod.veg16.boom, mod.veg17.boom, mod.veg18.boom, mod.veg19.boom)

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
aic.boom$model <- c("null", "grass*evi+elevation+fire", "grass+evi+elevation+fire", "grass*evi+fire", "grass*evi+elevation", "grass+evi+fire", "grass+evi+elevation", "grass*evi", "grass+evi", "grass+elevation+fire", "evi+elevation+fire", "grass+elevation", "grass+fire", "elevation+fire", "evi+fire", "elevation+evi", "grass", "evi", "elevation", "fire")
aic.boom$id <- names(mods.boom)
aic.boom <- aic.boom %>% 
  arrange(delta)

mods.call <- list(mod.null.call, mod.veg1.call, mod.veg2.call, mod.veg3.call, mod.veg4.call, mod.veg5.call, mod.veg6.call, mod.veg7.call, mod.veg8.call, mod.veg9.call, mod.veg10.call, mod.veg11.call, mod.veg12.call, mod.veg13.call, mod.veg14.call, mod.veg15.call, mod.veg16.call, mod.veg17.call, mod.veg18.call, mod.veg19.call)

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
aic.call$model <- c("null", "grass*evi+elevation+fire", "grass+evi+elevation+fire", "grass*evi+fire", "grass*evi+elevation", "grass+evi+fire", "grass+evi+elevation", "grass*evi", "grass+evi", "grass+elevation+fire", "evi+elevation+fire", "grass+elevation", "grass+fire", "elevation+fire", "evi+fire", "elevation+evi", "grass", "evi", "elevation", "fire")
aic.call <- aic.call %>% 
  arrange(delta)

mods.extra <- list(mod.null.extra, mod.veg1.extra, mod.veg2.extra, mod.veg3.extra, mod.veg4.extra, mod.veg5.extra, mod.veg6.extra, mod.veg7.extra, mod.veg8.extra, mod.veg9.extra, mod.veg10.extra, mod.veg11.extra, mod.veg12.extra, mod.veg13.extra, mod.veg14.extra, mod.veg15.extra, mod.veg16.extra, mod.veg17.extra, mod.veg18.extra, mod.veg19.extra)

aic.extra <- data.frame(df=sapply(mods.extra, function(z) length(coef(z))),
                       AIC=sapply(mods.extra, AIC),
                       loglik = sapply(mods.extra, function(z) logLik(z))) %>% 
  mutate(AICc = AIC + (2*df^2+2*df) / (nobs(mod.null.extra)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight) %>% 
  mutate(response = "extra")
aic.extra$model <- c("null", "grass*cover+elevation+fire", "grass+cover+elevation+fire", "grass*cover+fire", "grass*cover+elevation", "grass+cover+fire", "grass+cover+elevation", "grass*cover", "grass+cover", "grass+elevation+fire", "cover+elevation+fire", "grass+elevation", "grass+fire", "elevation+fire", "cover+fire", "elevation+cover", "grass", "cover", "elevation", "fire")
aic.extra <- aic.extra %>% 
  arrange(delta)

aic <- rbind(aic.boom, aic.call, aic.extra) %>% 
  arrange(response, delta) %>% 
  dplyr::select(model, df, loglik, AICc, delta, weight, response)
View(aic)

write.csv(aic, "OccupancyModelSelection.csv", row.names = FALSE)

summary(mod.veg11.boom)
best.boom <- mod.veg11.boom
saveRDS(best.boom, "OccupancyModel_Boom.rds")

summary(mod.veg11.call)
best.call <- mod.veg11.call
saveRDS(best.call, "OccupancyModel_Call.rds")

summary(mod.null.extra)
best.extra <- mod.null.extra
saveRDS(best.call, "OccupancyModel_Extra.rds")

#9. Model predictions----
#Newdata
newdat <- data.frame(expand_grid(grass = seq(0, 1, 0.01),
                                 elevation = seq(round(min(dat$Elevation)), round(max(dat$Elevation)), 1))) %>% 
  mutate(grass.s = (grass - mean(dat$grass.300))/sd(dat$grass.300),
         elevation.s = (elevation - mean(dat$Elevation))/sd(dat$Elevation))

#Predict
Xnewdat.boom <- model.matrix(~grass.s + elevation.s, newdat)
newdat$delta.boom <- plogis(drop(Xnewdat.boom %*% best.boom$coef[1:3]))
newdat$delta.boom.sd <- plogis(drop(Xnewdat.boom%*% best.boom$summary[1:3,2]))
newdat$delta.boom.high <- newdat$delta.boom + newdat$delta.boom.sd/sqrt(length(X.boom))*1.96
newdat$delta.boom.low <- newdat$delta.boom - newdat$delta.boom.sd/sqrt(length(X.boom))*1.96

Xnewdat.call <- model.matrix(~grass.s + elevation.s, newdat)
newdat$delta.call <- plogis(drop(Xnewdat.call %*% best.call$coef[1:3]))
newdat$delta.call.sd <- plogis(drop(Xnewdat.call%*% best.call$summary[1:3,2]))
newdat$delta.call.high <- newdat$delta.call + newdat$delta.call.sd/sqrt(length(X.call))*1.96
newdat$delta.call.low <- newdat$delta.call - newdat$delta.call.sd/sqrt(length(X.call))*1.96

write.csv(newdat, "OccupancyModelPredictions.csv", row.names = FALSE)

#Plot
newdat.grass <- newdat %>% 
  dplyr::filter(elevation==round(mean(dat$Elevation), 0))

newdat.elev <- newdat %>% 
  dplyr::filter(grass==round(mean(dat$grass.300), 2))

plot.boom.grass <- ggplot(newdat.grass) +
  geom_ribbon(aes(x=grass, ymin=delta.boom.low, ymax = delta.boom.high), alpha = 0.5) +
  geom_line(aes(x=grass, y=delta.boom), show.legend = FALSE)

plot.call.grass <- ggplot(newdat.grass) +
  geom_ribbon(aes(x=grass, ymin=delta.call.low, ymax = delta.call.high), alpha = 0.5) +
  geom_line(aes(x=grass, y=delta.call))

plot.boom.elev <- ggplot(newdat.elev) +
  geom_ribbon(aes(x=elevation, ymin=delta.boom.low, ymax = delta.boom.high), alpha = 0.5) +
  geom_line(aes(x=elevation, y=delta.boom))

plot.call.elev <- ggplot(newdat.elev) +
  geom_ribbon(aes(x=elevation, ymin=delta.call.low, ymax = delta.call.high), alpha = 0.5) +
  geom_line(aes(x=elevation, y=delta.call))

grid.arrange(plot.boom.grass, plot.call.grass,
             plot.boom.elev, plot.call.elev,
             nrow=2, ncol=2)

#10. Estimates----
best.boom$coef
best.call$coef
best.extra$coef

#10a. Suitability----
delta.boom <- data.frame(delta = plogis(drop(Xveg11.boom %*% best.boom$coef[1:3])),
                         sa=station.boom$station)
summary(delta.boom$delta)

delta.call <- data.frame(delta = plogis(drop(Xveg11.call %*% best.call$coef[1:3])),
                         sa=station.call$station)
summary(delta.call$delta)

#10b. Occupied----
summary(1-exp(-lambda$lambda[1]))
summary(1-exp(-lambda$lambda[2]))

#10c. Active----
summary(as.numeric(p.boom))
summary(as.numeric(p.call))

#10d. Present (closure)----
plogis(best.boom$coef[4])
plogis(best.boom$summary[4,2])

plogis(best.call$coef[4])
plogis(best.call$summary[4,2])

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
  

save.image("DEOccupancyModels.Rdata")
#load("DEOccupancyModels.Rdata")

