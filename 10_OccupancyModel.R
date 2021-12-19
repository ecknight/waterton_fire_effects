#title: ZIP abundance model fitting for CONI density modelling
#author: Elly C. Knight, adapted from Peter Solymos
#created: November 8, 2021

library(tidyverse)
library(lubridate)
library(tidylog)
library(AICcmodavg)
library(sf)
library(dggridR)
library(data.table)

options(scipen=99999)

source("src/mvocc.R")

#1. Wrangle----
dat.raw <- read.csv("SurveyDataWithCovs.csv") %>% 
  mutate(ID=paste0(survey,"-",station,"-",year),
         boom=ifelse(detection==2, 1, 0),
         call=ifelse(detection>0, 1, 0),
         DateTime = ymd_hms(DateTime),
         doy = yday(DateTime),
         FireTimeBin = ifelse(FireTime <=4, 1, 0)) %>% 
  arrange(ID, DateTime) %>% 
  group_by(ID) %>% 
  mutate(n=row_number()) %>% 
  ungroup() %>% 
  mutate(p.boom = p.boom/(max(p.boom)),
         p.call = p.call/max(p.call))

dat.occ <- dat.raw %>% 
  group_by(ID, survey) %>% 
  summarize(occ.boom = ifelse(sum(boom) > 0, 1, 0),
            occ.call = ifelse(sum(call) > 0, 1, 0)) %>% 
  ungroup()

table(dat.occ$occ.boom, dat.occ$occ.call, dat.occ$survey)

dat <- dat.raw %>% 
  left_join(dat.occ)

#2. Spatial separation of points----
dat.dist <- dat %>% 
  dplyr::select(station, X, Y) %>% 
  unique() %>% 
  st_as_sf(coords=c("X", "Y"), crs=26912) %>% 
  st_distance() %>% 
  data.frame() %>% 
  pivot_longer(cols=X1:X269, names_to="station", values_to="distance") %>% 
  mutate(distance = as.numeric(distance)) %>% 
  dplyr::filter(distance > 0, distance < 1000)

hist(dat.dist$distance)
#Ok yeah should probably look at spatially thinning - might solve the matrix issue too

#3. Set up sampling grid----
grid <- dgconstruct(area=0.5, metric=TRUE)

dat.station <- dat %>% 
  dplyr::select(station, X, Y) %>% 
  unique() %>% 
  st_as_sf(coords=c("X", "Y"), crs=26912) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(Latitude = Y, Longitude = X) %>% 
  cbind(dat %>% 
          dplyr::select(station, X, Y) %>% 
          unique())

dat.grid <- dat.station %>% 
  mutate(cell = dgGEO_to_SEQNUM(grid, Longitude, Latitude)$seqnum) %>% 
  arrange(cell) %>% 
  full_join(dat)

#4. Set up loop for habitat model----
boot <- 100
set.seed(999)
method <- "Nelder-Mead" # fail fast
#method <- "SANN"
#method <- "DE" # slow and sure

aic.list <- list()
for(i in 1:boot){
  
  dat.i <- dat.grid %>% 
    dplyr::select(station, cell) %>% 
    unique() %>% 
    group_by(cell) %>% 
#    sample_n(1) %>% 
    ungroup() %>% 
    left_join(dat.grid) 
  
  #5. Format for occupancy----
  lambda <- read.csv("LambdaEstimates.csv")
  
  #Boom
  dat.boom <- dat.i %>% 
    dplyr::filter(survey=="ARU")
  
  station.boom <- dat.boom %>% 
    dplyr::select(station, year, fire, ID, X, Y, Elevation, FireTime, grass.300, sand.300, cover.300, Elevation) %>% 
    unique() %>% 
    dplyr::filter(!(ID=="ARU-WLNP-1-1-2018" & cover.300>1.1290320))
  #figure out what's actually going on with this ID
  
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
  dat.call <- dat.i
  
  station.call <- dat.call %>% 
    dplyr::select(station, year, fire, ID, X, Y, Elevation, FireTime, grass.300, sand.300, trails.300, cover.300, Elevation) %>% 
    unique() %>% 
    dplyr::filter(!(ID=="ARU-WLNP-1-1-2018" & cover.300>1.1290320))
  
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
  
  #6. Model----
  #model matrices
  X.boom <- matrix(1, M.boom, 1)
  Xveg1.boom <- model.matrix(~Elevation*cover.300 + grass.300, station.boom)
  Xveg2.boom <- model.matrix(~Elevation*cover.300, station.boom)
  Xveg3.boom <- model.matrix(~Elevation + cover.300 + grass.300, station.boom)
  Xveg4.boom <- model.matrix(~grass.300 + cover.300, station.boom)
  Xveg5.boom <- model.matrix(~grass.300 + Elevation, station.boom)
  Xveg6.boom <- model.matrix(~cover.300 + Elevation, station.boom)
  Xveg7.boom <- model.matrix(~grass.300, station.boom)
  Xveg8.boom <- model.matrix(~cover.300, station.boom)
  Xveg9.boom <- model.matrix(~Elevation, station.boom)
  
  X.call <- matrix(1, M.call, 1)
  Xveg1.call <- model.matrix(~Elevation*cover.300 + grass.300, station.call)
  Xveg2.call <- model.matrix(~Elevation*cover.300, station.call)
  Xveg3.call <- model.matrix(~Elevation + cover.300 + grass.300, station.call)
  Xveg4.call <- model.matrix(~grass.300 + cover.300, station.call)
  Xveg5.call <- model.matrix(~grass.300 + Elevation, station.call)
  Xveg6.call <- model.matrix(~cover.300 + Elevation, station.call)
  Xveg7.call <- model.matrix(~grass.300, station.call)
  Xveg8.call <- model.matrix(~cover.300, station.call)
  Xveg9.call <- model.matrix(~Elevation, station.call)

  #Boom
  mod.null.boom <- try(mvocc(y.boom, X.boom, Z.boom, p.boom, lamvec.boom, method=method))
  mod.null.call <- try(mvocc(y.call, X.call, Z.call, p.call, lamvec.call, method=method))
  
  
  if(class(mod.null.boom)!="mvocc"){
    next
  }
  else{
    mod.veg1.boom <- mvocc(y.boom, Xveg1.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg2.boom <- mvocc(y.boom, Xveg2.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg3.boom <- mvocc(y.boom, Xveg3.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg4.boom <- mvocc(y.boom, Xveg4.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg5.boom <- mvocc(y.boom, Xveg5.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg6.boom <- mvocc(y.boom, Xveg6.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg7.boom <- mvocc(y.boom, Xveg7.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg8.boom <- mvocc(y.boom, Xveg8.boom, Z.boom, p.boom, lamvec.boom, method=method)
    mod.veg9.boom <- mvocc(y.boom, Xveg9.boom, Z.boom, p.boom, lamvec.boom, method=method)
    
    mod.veg1.call <- mvocc(y.call, Xveg1.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg2.call <- mvocc(y.call, Xveg2.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg3.call <- mvocc(y.call, Xveg3.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg4.call <- mvocc(y.call, Xveg4.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg5.call <- mvocc(y.call, Xveg5.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg6.call <- mvocc(y.call, Xveg6.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg7.call <- mvocc(y.call, Xveg7.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg8.call <- mvocc(y.call, Xveg8.call, Z.call, p.call, lamvec.call, method=method)
    mod.veg9.call <- mvocc(y.call, Xveg9.call, Z.call, p.call, lamvec.call, method=method)
    
    aic.boom <- AIC(mod.null.boom, mod.veg1.boom, mod.veg2.boom, mod.veg3.boom, mod.veg4.boom, mod.veg5.boom, mod.veg6.boom, mod.veg7.boom, mod.veg8.boom, mod.veg9.boom) %>% 
      arrange(AIC)
    aic.boom$model <- row.names(aic.boom)
    aic.boom$response <- "boom"
    
    aic.call <- AIC(mod.null.call, mod.veg1.call, mod.veg2.call, mod.veg3.call, mod.veg4.call, mod.veg5.call, mod.veg6.call, mod.veg7.call, mod.veg8.call, mod.veg9.call) %>% 
      arrange(AIC)
    aic.call$model <- row.names(aic.call)
    aic.call$response <- "call"
    
    aic.list[[i]] <- rbind(aic.boom, aic.call) %>% 
#      mutate(boot = i)
    
    aic.list[[i]] <- aic.boom %>% 
      mutate(boot = i)
    
#    print(paste0("Finished bootstrap ", i, " of ", boot, " - ", head(aic.boom$model), " & ", head(aic.call$model)))
    print(paste0("Finished bootstrap ", i, " of ", boot, " - ", head(aic.boom$model, 1)))
    
  }

  
}  

aic <- rbindlist(aic.list) %>% 
  group_by(model, response) %>% 
  summarize(mean = mean(AIC), 
            sd = sd(AIC)) %>% 
  ungroup() %>% 
  arrange(response, mean)
View(aic)
  




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
