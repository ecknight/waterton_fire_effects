#title: Covariate selection via boosted regression tree for CONI density modelling
#author: Elly C. Knight
#created: November 7, 2021

library(tidyverse)
library(lubridate)
library(dismo)

#idea: only use ARU data and summarize at station level

#1. Wrangle----
dat <- read.csv("SurveyDataWithCovs.csv") %>% 
  mutate(ID=paste0(survey,"-",station,"-",year),
         boom=ifelse(detection==2, 1, 0),
         call=ifelse(detection %in% c(1,2,3), 1, 0),
         extra=ifelse(detection==3, 1, 0),
         DateTime = ymd_hms(DateTime),
         doy = yday(DateTime)) %>% 
  dplyr::filter(survey=="ARU") %>% 
  group_by(ID, Elevation, cover.300, develop.300, grass.300, pine.300, water.300, wetland.300, evi) %>% 
  summarize(boom = ifelse(sum(boom) > 0, 1, 0),
            call= ifelse(sum(call) > 0, 1, 0),
            extra = ifelse(sum(extra) > 0, 1, 0),
            p.boom = mean(p.boom),
            p.call = mean(p.call),
            p.extra = mean(p.extra)) %>% 
  ungroup() %>% 
  dplyr::select(boom, p.boom, call, p.call, extra, p.extra, Elevation, cover.300, develop.300, grass.300,  pine.300, water.300, wetland.300, evi) %>% 
  data.frame() %>% 
  dplyr::filter(Elevation > 0)

#3. Model boom----
set.seed(1234)
boom.gbm <- dismo::gbm.step(data=dat, 
                     gbm.x=7:14,
                     gbm.y=1,
                     offset=dat$p.boom,
                     family="bernoulli",
                     tree.complexity = 3,
                     learning.rate = 0.001,
                     bag.fraction = 0.75,
                     max.trees=10000) 

gbm.plot(boom.gbm)

boom.int <- gbm.interactions(boom.gbm)

#4. Model call----
set.seed(1234)
call.gbm <- dismo::gbm.step(data=dat, 
                            gbm.x=7:14,
                            gbm.y=3,
                            offset=dat$p.call,
                            family="bernoulli",
                            tree.complexity = 3,
                            learning.rate = 0.001,
                            bag.fraction = 0.75,
                            max.trees=10000) 

gbm.plot(call.gbm)

call.int <- gbm.interactions(call.gbm)

#5. Model extraterritorial----
set.seed(1234)
extra.gbm <- dismo::gbm.step(data=dat, 
                            gbm.x=7:14,
                            gbm.y=5,
                            offset=dat$p.extra,
                            family="bernoulli",
                            tree.complexity = 3,
                            learning.rate = 0.001,
                            bag.fraction = 0.75,
                            max.trees=10000) 

gbm.plot(extra.gbm)

extra.int <- gbm.interactions(extra.gbm)

#5. Look at results----
summary <- data.frame(summary(boom.gbm)) %>% 
  mutate(response="boom") %>% 
  rbind(data.frame(summary(call.gbm)) %>% 
          mutate(response="call")) %>% 
  rbind(data.frame(summary(extra.gbm)) %>% 
          mutate(response="extra"))

View(summary)

write.csv(summary, "BRTCovariateResults.csv", row.names = FALSE)

int <- data.frame(boom.int$rank.list) %>% 
  mutate(response="boom") %>% 
  rbind(data.frame(call.int$rank.list) %>% 
          mutate(response="call")) %>% 
  rbind(data.frame(extra.int$rank.list) %>% 
          mutate(response="extra"))
int

write.csv(int, "BRTInteractionResults.csv", row.names = FALSE)
