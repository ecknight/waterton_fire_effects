#title: Covariate selection via boosted regression tree for CONI density modelling
#author: Elly C. Knight
#created: November 7, 2021

library(tidyverse)
library(lubridate)
library(dismo)

#1. Wrangle----
dat <- read.csv("SurveyDataWithOffsets.csv")  %>% 
  mutate(ID=paste0(survey,"-",station,"-",year),
         boom=ifelse(detection==2, 1, 0),
         call=ifelse(detection>0, 1, 0),
         DateTime = ymd_hms(DateTime),
         doy = yday(DateTime)) %>% 
  arrange(ID, DateTime) %>% 
  group_by(ID) %>% 
  mutate(n=row_number()) %>% 
  ungroup() %>% 
  dplyr::filter(doy < 225)

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
         FireSeverity = ifelse(year > 2017, FireSeverity, NA)) %>% 
  dplyr::filter(!is.na(pine.300))

write.csv(dat.cov, "SurveyDataWithCovs.csv", row.names = FALSE)

dat.gbm <- dat.cov %>% 
  dplyr::select(boom, p.boom, call, p.call, Elevation, cover.300, develop.300, grass.300, trails.300, pine.300, sand.300, water.300, wetland.300) %>% 
  data.frame()

#3. Model boom----

boom.gbm <- dismo::gbm.step(data=dat.gbm, 
                     gbm.x=5:13,
                     gbm.y=1,
                     offset=dat.gbm$p.boom,
                     family="bernoulli",
                     tree.complexity = 3,
                     learning.rate = 0.001,
                     bag.fraction = 0.75,
                     max.trees=10000) 

boom.int <- gbm.interactions(boom.gbm)

#4. Model call----
call.gbm <- dismo::gbm.step(data=dat.gbm, 
                            gbm.x=5:13,
                            gbm.y=3,
                            offset=dat.gbm$p.call,
                            family="bernoulli",
                            tree.complexity = 3,
                            learning.rate = 0.005,
                            bag.fraction = 0.75,
                            max.trees=10000) 

call.int <- gbm.interactions(call.gbm)

#5. Look at results----
summary(boom.gbm)
summary(call.gbm)

boom.int$rank.list
call.int$rank.list
