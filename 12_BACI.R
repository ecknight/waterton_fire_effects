library(tidyverse)
library(lubridate)
library(lme4)
library(MuMIn)
library(AICcmodavg)
library(suncalc)
library(sf)

options(scipen=99999)

#1. Wrangling----
dat <- read.csv("SurveyDataWithOffsets&Covariates.csv") %>% 
  mutate(Kenow = ifelse(FireHistory==2017, "impact", "control"),
         Kenow = ifelse(is.na(Kenow), "control", Kenow),
         date = date(DateTime))

dat.wgs <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=26912) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lon=X, lat=Y) %>% 
  cbind(dat) %>% 
  dplyr::select(date, lat, lon)

dat.sun <- getSunlightTimes(data=dat.wgs, keep="sunset") %>% 
  dplyr::select(sunset)
head(dat.sun)

dat.human <- dat %>% 
  cbind(dat.sun) %>% 
  dplyr::filter(survey=="human") %>% 
  mutate(DateTime = ymd_hms(DateTime), 
         doy = yday(DateTime),
         suntime = as.numeric(DateTime - sunset + 6)) %>% 
  dplyr::filter(doy < 210)

hist(dat.human$suntime)
hist(dat.human$doy)

#2. Data availability----
table(dat$fire, dat$Kenow, dat$survey)
table(dat.human$station, dat.human$fire, dat.human$Kenow)

#3. Ratios----
dat.human.sum <- dat.human %>% 
  group_by(fire, Kenow) %>% 
  summarize(n = n(),
            det = sum(detection)) %>% 
  ungroup() %>% 
  mutate(ratio = det/n)
dat.human.sum

#4. Temporal stuff---
ggplot(dat.human) +
  geom_jitter(aes(x=doy, y=detection)) +
  geom_smooth(aes(x=doy, y=detection))

ggplot(dat.human) +
  geom_jitter(aes(x=suntime, y=detection)) +
  geom_smooth(aes(x=suntime, y=detection))

#4. Model----
mod.int <- glmer(detection ~ fire*Kenow + poly(suntime, 2) + (1|station),
                 data=dat.human, family="binomial", na.action = "na.fail")
mod.add <- glmer(detection ~ fire + Kenow +  poly(suntime, 2) + (1|station),
                 data=dat.human, family="binomial", na.action = "na.fail")
mod.fire <- glmer(detection ~ fire +  poly(suntime, 2) + (1|station),
                  data=dat.human, family="binomial", na.action = "na.fail")
mod.kenow <- glmer(detection ~ Kenow +  poly(suntime, 2) + (1|station),
                   data=dat.human, family="binomial", na.action = "na.fail")
mod.null <- glmer(detection ~ 1  +  poly(suntime, 2) + (1|station),
                  data=dat.human, family="binomial", na.action = "na.fail")

mods <- list("mod.int" = mod.int, "mod.add" = mod.add, "mod.fire" = mod.fire, "mod.kenow" = mod.kenow, "mod.null" = mod.null)
aictab(mods)

summary(mod.null)
