#title: Covariate selection via boosted regression tree for CONI density modelling
#author: Elly C. Knight
#created: November 7, 2021

library(tidyverse)
library(lubridate)
library(dismo)

#1. Wrangle----
dat <- read.csv("SurveyDataWithCovs.csv") %>% 
  mutate(ID=paste0(survey,"-",station,"-",year),
         boom=ifelse(detection==2, 1, 0),
         call=ifelse(detection>0, 1, 0),
         DateTime = ymd_hms(DateTime),
         doy = yday(DateTime))

dat.gbm <- dat %>% 
  dplyr::select(boom, p.boom, call, p.call, Elevation, cover.300, develop.300, grass.300, trails.300, pine.300, sand.300, water.300, wetland.300) %>% 
  data.frame()

#3. Model boom----
set.seed(1234)
boom.gbm <- dismo::gbm.step(data=dat.gbm, 
                     gbm.x=5:13,
                     gbm.y=1,
                     offset=dat.gbm$p.boom,
                     family="bernoulli",
                     tree.complexity = 3,
                     learning.rate = 0.005,
                     bag.fraction = 0.75,
                     max.trees=10000) 

boom.int <- gbm.interactions(boom.gbm)

#4. Model call----
set.seed(1234)
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
summary <- data.frame(summary(boom.gbm)) %>% 
  mutate(response="boom") %>% 
  rbind(data.frame(summary(call.gbm)) %>% 
          mutate(response="call"))

write.csv(summary, "BRTCovariateResults.csv", row.names = FALSE)

int <- data.frame(boom.int$rank.list) %>% 
  mutate(response="boom") %>% 
  rbind(data.frame(call.int$rank.list) %>% 
          mutate(response="call"))

write.csv(summary, "BRTInteractionResults.csv", row.names = FALSE)