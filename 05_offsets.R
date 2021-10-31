#title: Planning for 2021 ARU deployment to fill sampling gaps for CONI
#author: Elly C. Knight, adapted from Peter Solymos
#created: October 31, 2021

library(mefa4)
library(maptools)
library(survival)
library(tidyverse)
library(readr)
library(fs)
library(lubridate)

#This step uses just the new recognizer data to calculate availability offsets for date and time

#1. Read in data----
df.raw <- read.csv("Data/ValidatedRecognizerResults.csv") %>% 
  mutate(timeofdetection=as.numeric(str_sub(TimeOffset, 5, 5))*60+as.numeric(str_sub(TimeOffset, 7, 8))) %>% 
  mutate(DateTime = as_datetime(DateTime))

#2. Add files with no detections----
rec <- read.csv("FinalFileList.csv") %>% 
  separate(path, 
           into = c("Blank", "Volumes", "Drive", "Park","Batch", "Filename"),
           sep = "/",
           remove = FALSE) %>% 
  dplyr::filter(!is.na(Filename),
                !Filename %in% df$Station) %>% 
  separate(Filename, 
           into = c("Station", "Date", "Time", "Filetype"),
           sep = "[^A-Za-z0-9 \\- :]",
           remove = FALSE) %>% 
  mutate(detection="0",
         TimeOffset=NA,
         timeofdetection=NA,
         Date = ymd(Date),
         DateTime = as_datetime(paste0(Date, Time)), .after = Time) %>% 
  dplyr::select(Filename, Station, Date, Time, DateTime, TimeOffset, detection, timeofdetection)

df <- rbind(df.raw, rec)

#3. ID locations with detections----
site.boom <- df %>% 
  filter(detection=="B") %>% 
  dplyr::select(Station) %>% 
  unique()

df.boom <- df %>% 
  filter(Station %in% site.boom$Station) %>% 
  mutate(detection = ifelse(detection=="B", 1, 0))

site.call <- df %>% 
  filter(detection %in% c("1", "B")) %>% 
  dplyr::select(Station) %>% 
  unique()

#4. Subset to time to first detection and treat nondetections as censored events----
#also make 0 seconds 1 seconds
df.boom.1st <- df %>% 
  filter(Station %in% site.boom$Station) %>% 
  mutate(detection = ifelse(detection=="B", 1, 0)) %>% 
  arrange(Filename, timeofdetection) %>% 
  group_by(Filename) %>% 
  mutate(order=row_number()) %>% 
  ungroup() %>% 
  dplyr::filter(order==1) %>% 
  dplyr::select(-order) %>% 
  mutate(time = case_when(is.na(timeofdetection) ~ 600,
                                     timeofdetection==0 ~ 1,
                                     !is.na(timeofdetection) ~ timeofdetection))

df.call.1st <- df %>% 
  filter(Station %in% site.call$Station) %>% 
  mutate(detection = as.numeric(ifelse(detection=="B", 1, detection)))  %>% 
  arrange(Filename, timeofdetection) %>% 
  group_by(Filename) %>% 
  mutate(order=row_number()) %>% 
  ungroup() %>% 
  dplyr::filter(order==1) %>% 
  dplyr::select(-order) %>% 
  mutate(time = case_when(is.na(timeofdetection) ~ 600,
                                     timeofdetection==0 ~ 1,
                                     !is.na(timeofdetection) ~ timeofdetection))

#5. Add additional fields----
df.boom <- df.boom.1st %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         tod = start/24,
         sin = sin(tod*2*pi),
         cos = cos(tod*2*pi),
         sv = Surv(time/60, detection))

df.call <- df.call.1st %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         tod = start/24,
         sin = sin(tod*2*pi),
         cos = cos(tod*2*pi),
         sv = Surv(time/60, detection))

#6. Check for file duplicates----
stopifnot(!any(duplicated(df.boom$Filename)))
stopifnot(!any(duplicated(df.call$Filename)))

#7. Mean of survival time in minutes----
mean(df.boom$sv)
hist(df.boom$time/60)
mean(df.call$sv)
hist(df.call$time/60)

#8. Fit a series of models----
mods.boom <- list(
  m0 = survreg(sv ~ 1, df.boom, dist="exponential"),
  m1 = survreg(sv ~ jday, df.boom, dist="exponential"),
  m2 = survreg(sv ~ jday + jday2, df.boom, dist="exponential"),
  m3 = survreg(sv ~ cos, df.boom, dist="exponential"),
  m4 = survreg(sv ~ cos + sin, df.boom, dist="exponential"),
  m5 = survreg(sv ~ jday + cos, df.boom, dist="exponential"),
  m6 = survreg(sv ~ jday + cos + sin, df.boom, dist="exponential"),
  m7 = survreg(sv ~ jday + jday2 + cos, df.boom, dist="exponential"),
  m8 = survreg(sv ~ jday + jday2 + cos + sin, df.boom, dist="exponential")
)


mods.call <- list(
  m0 = survreg(sv ~ 1, df.call, dist="exponential"),
  m1 = survreg(sv ~ jday, df.call, dist="exponential"),
  m2 = survreg(sv ~ jday + jday2, df.call, dist="exponential"),
  m3 = survreg(sv ~ cos, df.call, dist="exponential"),
  m4 = survreg(sv ~ cos + sin, df.call, dist="exponential"),
  m5 = survreg(sv ~ jday + cos, df.call, dist="exponential"),
  m6 = survreg(sv ~ jday + cos + sin, df.call, dist="exponential"),
  m7 = survreg(sv ~ jday + jday2 + cos, df.call, dist="exponential"),
  m8 = survreg(sv ~ jday + jday2 + cos + sin, df.call, dist="exponential")
)

#9. Compare with AIC----
aic.boom <- data.frame(df=sapply(mods.boom, function(z) length(coef(z))),
                  AIC=sapply(mods.boom, AIC),
                  loglik = sapply(mods.boom, function(z) logLik(z))) %>%
  mutate(AICc = AIC + (2*df^2+2*df) / (length(df.boom)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight) %>%
  arrange(delta)
aic.boom


aic.call <- data.frame(df=sapply(mods.call, function(z) length(coef(z))),
                       AIC=sapply(mods.call, AIC),
                       loglik = sapply(mods.call, function(z) logLik(z))) %>%
  mutate(AICc = AIC + (2*df^2+2*df) / (length(df.call)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight) %>%
  arrange(delta)
aic.call

#Seems like cos + sin fits best for both