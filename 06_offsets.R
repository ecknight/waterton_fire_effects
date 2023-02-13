#title: Availability for detection offsets for CONI density modelling
#author: Elly C. Knight, adapted from Peter Solymos
#created: October 31, 2021

library(mefa4)
library(maptools)
library(survival)
library(tidyverse)
#library(tidylog)
library(readr)
library(fs)
library(lubridate)
library(suncalc)

options(scipen=99999)

#This step uses just the new recognizer data to model probability of availability for date and time and then uses it to predict offsets for human and recognizer data that will be used in the models

#1. Read in recognizer data----
df.raw <- read.csv("Data/ValidatedRecognizerResults.csv") %>% 
  mutate(timeofdetection=as.numeric(str_sub(TimeOffset, 5, 5))*60+as.numeric(str_sub(TimeOffset, 7, 8))) %>%
  rename(path = Filename) %>% 
  mutate(DateTime = as_datetime(DateTime),
         year = year(DateTime),
         path = gsub("\\", "/", path, fixed=TRUE),
         path = gsub("/Volumes/", "", path, fixed=TRUE)) %>% 
  separate(path, 
           into = c("Drive", "Park", "Batch", "Filename"),
           sep = "/",
           remove = FALSE) %>% 
  mutate(Date = as.character(Date)) %>% 
  dplyr::select(Filename, Station, year, Date, Time, DateTime, TimeOffset, detection, timeofdetection) %>% 
  dplyr::filter(detection!="0")

#2. Add files with no detections----
rec <- read.csv("FinalFileList.csv") %>% 
  separate(path, 
           into = c("Blank", "Volumes", "Drive", "Park","Batch", "Filename"),
           sep = "/",
           remove = FALSE) %>% 
  dplyr::filter(!is.na(Filename)) %>% 
  separate(Filename, 
           into = c("Station", "Date", "Time", "Filetype"),
           sep = "[^A-Za-z0-9 \\- :]",
           remove = FALSE) %>% 
  mutate(detection="0",
         TimeOffset=NA,
         timeofdetection=NA,
         Date = ymd(Date),
         Time = as.character(Time),
         DateTime = ymd_hms(paste0(Date, Time)),
         year = year(DateTime)) %>% 
  dplyr::select(Filename, Station, year, Date, Time, DateTime, TimeOffset, detection, timeofdetection) %>% 
  dplyr::filter(!Filename %in% df.raw$Filename,
                Station!="WLNP-022-1", #remove files for duplicate station ID
                !Station %in% c("WLNP-023-04", "WLNP-13-3", "WLNP-16-1", "WLNP-8-3", "WLNP-83")) #remove files that don't have 30 recordings
#Remember: Missing pieces are from folders, not recordings

df <- rbind(df.raw, rec)

tz(df$DateTime) <- "Canada/Mountain"

hist(hour(df$DateTime))

#3. Add time since sunset-----
#Use mean location, doesn't really matter
locs <- read.csv("ExistingSamplingLocations.csv") %>% 
  dplyr::select(X, Y) %>% 
  summarize(X=mean(X),
            Y=mean(Y))

df.locs <- df %>% 
  mutate(lon=locs$X,
         lat=locs$Y) %>% 
  mutate(date=ymd(str_sub(as.character(DateTime), 1, 10)))

df.sun <- getSunlightTimes(data=df.locs, tz="Canada/Mountain", keep="sunset") %>% 
  cbind(df) %>% 
  mutate(tsss = as.numeric(difftime(DateTime, sunset, units="hours")),
         tsss = ifelse(tsss < -12, tsss+24, tsss)) 

ggplot(df.sun) +
  geom_histogram(aes(x=tsss))

write.csv(df.sun, "ZeroFilledRecognizerResults.csv", row.names = FALSE)

#4. ID locations with detections----

df.sun <- read.csv("ZeroFilledRecognizerResults.csv")

site.boom <- df.sun %>% 
  filter(detection=="B") %>% 
  dplyr::select(Station) %>% 
  unique()

site.call <- df.sun %>% 
  filter(detection %in% c("1", "B")) %>% 
  dplyr::select(Station) %>% 
  unique()

site.extra <- site.call %>% 
  dplyr::filter(!Station %in% site.boom$Station)


#5. Subset to time to first detection and treat nondetections as censored events----
#also make 0 seconds 1 seconds
df.boom.1st <- df.sun %>% 
  filter(Station %in% site.boom$Station) %>% 
  mutate(detection = ifelse(detection=="B", 1, 0),
         timeofdetection = ifelse(detection==1, timeofdetection, 600)) %>% 
  arrange(Filename, timeofdetection) %>% 
  group_by(Filename) %>% 
  mutate(order=row_number()) %>% 
  ungroup() %>% 
  dplyr::filter(order==1) %>% 
  dplyr::select(-order) %>% 
  mutate(time = case_when(is.na(timeofdetection) ~ 600,
                                     timeofdetection==0 ~ 1,
                                     !is.na(timeofdetection) ~ timeofdetection))

df.call.1st <- df.sun %>% 
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

df.extra.1st <- df.sun %>% 
  filter(Station %in% site.extra$Station) %>% 
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

#6. Add additional fields----
df.boom <- df.boom.1st %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         tod = start/24,
         ts = (tsss +2)/8,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi),
         sv = Surv(time/60, detection))

ggplot(df.boom) +
  geom_point(aes(x=tsss, y=sin))

ggplot(df.boom) +
  geom_point(aes(x=tsss, y=cos))

ggplot(df.boom) +
  geom_point(aes(x=cos, y=cos + sin))

df.call <- df.call.1st %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         tod = start/24,
         ts = (tsss +2)/8,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi),
         sv = Surv(time/60, detection))

df.extra <- df.extra.1st %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         tod = start/24,
         ts = (tsss+2)/8,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi),
         sv = Surv(time/60, detection))

#7. Check for file duplicates----
stopifnot(!any(duplicated(df.boom$Filename)))
stopifnot(!any(duplicated(df.call$Filename)))
stopifnot(!any(duplicated(df.extra$Filename)))

#8. Mean of survival time in minutes----
mean(df.boom$sv)
hist(df.boom$time/60)
mean(df.call$sv)
hist(df.call$time/60)
mean(df.extra$sv)
hist(df.extra$time/60)

#9. Fit a series of models----
mods.boom <- list(
  m0 = survreg(sv ~ 1, df.boom, dist="exponential"),
  m1 = survreg(sv ~ jday, df.boom, dist="exponential"),
  m2 = survreg(sv ~ sin, df.boom, dist="exponential"),
  m3 = survreg(sv ~ cos, df.boom, dist="exponential"),
  m4 = survreg(sv ~ sin + cos, df.boom, dist="exponential"),
  m5 = survreg(sv ~ jday + sin, df.boom, dist="exponential"),
  m6 = survreg(sv ~ jday + cos, df.boom, dist="exponential"),
  m7 = survreg(sv ~ jday + sin + cos, df.boom, dist="exponential"))

mods.call <- list(
  m0 = survreg(sv ~ 1, df.call, dist="exponential"),
  m1 = survreg(sv ~ jday, df.call, dist="exponential"),
  m2 = survreg(sv ~ sin, df.call, dist="exponential"),
  m3 = survreg(sv ~ cos, df.call, dist="exponential"),
  m4 = survreg(sv ~ sin + cos, df.call, dist="exponential"),
  m5 = survreg(sv ~ jday + sin, df.call, dist="exponential"),
  m6 = survreg(sv ~ jday + cos, df.call, dist="exponential"),
  m7 = survreg(sv ~ jday + sin + cos, df.call, dist="exponential"))

mods.extra <- list(
  m0 = survreg(sv ~ 1, df.extra, dist="exponential"),
  m1 = survreg(sv ~ jday, df.extra, dist="exponential"),
  m2 = survreg(sv ~ sin, df.extra, dist="exponential"),
  m3 = survreg(sv ~ cos, df.extra, dist="exponential"),
  m4 = survreg(sv ~ sin + cos, df.extra, dist="exponential"),
  m5 = survreg(sv ~ jday + sin, df.extra, dist="exponential"),
  m6 = survreg(sv ~ jday + cos, df.extra, dist="exponential"),
  m7 = survreg(sv ~ jday + sin + cos, df.extra, dist="exponential"))

#10. Compare with AIC----
aic.boom <- data.frame(df=sapply(mods.boom, function(z) length(coef(z))),
                  AIC=sapply(mods.boom, AIC),
                  loglik = sapply(mods.boom, function(z) logLik(z))) %>%
  mutate(AICc = AIC + (2*df^2+2*df) / (length(df.boom)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight)
aic.boom$model <- c("1", "jday", "sin", "cos", "sin + cos", "jday + sin", "jday + cos", "jday + sin + cos")
aic.boom

aic.call <- data.frame(df=sapply(mods.call, function(z) length(coef(z))),
                       AIC=sapply(mods.call, AIC),
                       loglik = sapply(mods.call, function(z) logLik(z))) %>%
  mutate(AICc = AIC + (2*df^2+2*df) / (length(df.call)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight)
aic.call$model <- c("1", "jday", "sin", "cos", "sin + cos", "jday + sin", "jday + cos", "jday + sin + cos")
aic.call

aic.extra <- data.frame(df=sapply(mods.extra, function(z) length(coef(z))),
                       AIC=sapply(mods.extra, AIC),
                       loglik = sapply(mods.extra, function(z) logLik(z))) %>%
  mutate(AICc = AIC + (2*df^2+2*df) / (length(df.extra)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight)
aic.extramodel <- c("1", "jday", "sin", "cos", "sin + cos", "jday + sin", "jday + cos", "jday + sin + cos")
aic.extra

aic <- rbind(aic.boom %>% 
               mutate(response="boom"),
             aic.call %>% 
               mutate(response="call")) %>% 
  dplyr::select(model, df, loglik, AICc, delta, weight, response) %>% 
  mutate(loglik = round(loglik, 2),
         AICc = round(AICc, 2),
         delta = round(delta, 2),
         weight = round(weight, 2)) %>% 
  arrange(response, delta)

write.csv(aic, "OffsetAICresults.csv", row.names = FALSE)

#11. Summary of predictions----
#Boom
mb.boom <- mods.boom[[which.min(aic.boom$AICc)]]
summary(mb.boom)
#' Survival times
summary(predict(mb.boom))
#' Event rate per unit (1 min) time
summary(1/predict(mb.boom))
#' Probability of at least 1 event per 10 time units (10 mins)
summary(1-exp(-(1/predict(mb.boom))*10))

#Call
mb.call <- mods.call[[which.min(aic.call$AICc)]]
summary(mb.call)
#' Survival times
summary(predict(mb.call))
#' Event rate per unit (1 min) time
summary(1/predict(mb.call))
#' Probability of at least 1 event per 10 time units (10 mins)
summary(1-exp(-(1/predict(mb.call))*10))

#Extra
mb.extra <- mods.extra[[which.min(aic.extra$AICc)]]
summary(mb.extra)
#' Survival times
summary(predict(mb.extra))
#' Event rate per unit (1 min) time
summary(1/predict(mb.extra))
#' Probability of at least 1 event per 10 time units (10 mins)
summary(1-exp(-(1/predict(mb.extra))*10))

#12. Make predictions for figure----
#Boom
pr.boom <- data.frame(expand.grid(jday = seq(round(min(df.boom$jday), 3), round(max(df.boom$jday), 3), 0.001),
                      ts = seq(round(min(df.boom$ts), 2), round(max(df.boom$ts), 2), 0.01))) %>% 
  mutate(jday2 = jday^2,
         ts2 = ts^2,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi))
pr.boom$p.boom <- 1/predict(mb.boom, newdata=pr.boom)

#Call
pr.call <- data.frame(expand.grid(jday = seq(round(min(df.call$jday), 3), round(max(df.call$jday), 3), 0.001),
                                  ts = seq(round(min(df.call$ts), 2), round(max(df.call$ts), 2), 0.01))) %>% 
  mutate(jday2 = jday^2,
         ts2 = ts^2,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi))
pr.call$p.call <- 1/predict(mb.call, newdata=pr.call)

#Extra
pr.extra <- data.frame(expand.grid(jday = seq(round(min(df.extra$jday), 3), round(max(df.extra$jday), 3), 0.001),
                                  ts = seq(round(min(df.extra$ts), 2), round(max(df.extra$ts), 2), 0.01))) %>% 
  mutate(jday2 = jday^2,
         ts2 = ts^2,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi))
pr.extra$p.extra <- 1/predict(mb.extra, newdata=pr.extra)


preds <- full_join(pr.boom, pr.call) %>% 
  full_join(pr.extra) %>% 
  mutate(tsss = ts*8-2,
         doy = jday*365,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi))

write.csv(preds, "OffsetPredictionsForFigure.csv", row.names = FALSE)
write.csv(df.boom, "OffsetDataForFigure_Boom.csv", row.names = FALSE)
write.csv(df.call, "OffsetDataForFigure_Call.csv", row.names = FALSE)
write.csv(df.extra, "OffsetDataForFigure_Extra.csv", row.names = FALSE)

#13. Calculate offsets for recognizer data----
df.aru <- df.sun %>% 
  mutate(detection = ifelse(detection=="B", 2, as.numeric(detection))) %>% 
  group_by(Filename, Station, DateTime, tsss) %>% 
  summarize(detection=max(detection)) %>% 
  ungroup() %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         ts = (tsss +2)/8,
         ts2 = ts^2,
         sin = ts*2*pi,
         cos = ts*2*pi,
         survey="ARU",
         year = year(DateTime)) %>% 
  rename(station=Station)
df.aru$p.boom <- 1-exp(-(1/predict(mb.boom, newdata = df.aru))*10)
df.aru$p.call <- 1-exp(-(1/predict(mb.call, newdata = df.aru))*10)
df.aru$p.extra <- 1-exp(-(1/predict(mb.extra, newdata = df.aru))*10)

#14. Calculate offsets for human data----
df.hum.raw <- read.csv("Data/HumanData.csv") %>% 
  dplyr::filter(!is.na(Start.time)) %>% 
  mutate(Start.time=ymd_hm(paste0(Date, Start.time)),
         End.time=ymd_hm(paste0(Date, End.time)),
         Duration=End.time-Start.time) %>% 
  group_by(Year, Duration) %>% 
  mutate(Stops = n(),
         ID = row_number()) %>% 
  ungroup() %>% 
  mutate(Interval=Duration/Stops,
         DateTime = Start.time + Interval*ID,
         lon=locs$X,
         lat=locs$Y,
         date=ymd(str_sub(as.character(DateTime), 1, 10)))

tz(df.hum.raw$DateTime) <- "Canada/Mountain"

df.hum.sun <- getSunlightTimes(data=df.hum.raw, tz="Canada/Mountain", keep="sunset") %>% 
  cbind(df.hum.raw %>% 
          dplyr::select(-date, -lat, -lon)) %>% 
  mutate(tsss = as.numeric(difftime(DateTime, sunset, units="hours")),
         tsss = ifelse(tsss < -12, tsss+24, tsss)) 

ggplot(df.hum.sun) +
  geom_histogram(aes(x=tsss))

df.hum <- df.hum.sun %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         ts = (tsss +2)/8,
         ts2 = ts^2,
         sin = ts*2*pi,
         cos = ts*2*pi,
         survey="human") %>% 
  rename(station=Site, detection=Presence)

write.csv(df.hum, "HumanResults.csv", row.names=FALSE)

df.hum$p.boom <- NA
df.hum$p.call <- 1-exp(-(1/predict(mb.call, newdata = df.hum))*6)
df.hum$p.extra <- NA
         
#15. Correct for differences in recorder & recognizer EDR----
#From Yip et al. Fig 2A SM2 @ 4000 Hz
SM2 <- 0.789 

#Recall at score 60 from known distance clips
rec <-  0.5626959

df.hum$p.call <- df.hum$p.call*(1/SM2)*(1/rec)

summary(df.hum$p.boom)
summary(df.hum$p.call)

#15. Put together----
df.pred <- rbind(df.aru %>% 
                   dplyr::select(survey, station, DateTime, tsss, detection, p.boom, p.call, p.extra),
                 df.hum  %>% 
                   dplyr::select(survey, station, DateTime, tsss, detection, p.boom, p.call, p.extra)) %>%
  mutate(year=year(DateTime),
         doy=yday(DateTime),
         fire=ifelse(year <=2017, "before", "after")) %>% 
  dplyr::filter(doy < 225) %>% 
  data.frame() %>% 
  mutate(p.call = p.call/(max(p.call))) %>% 
  mutate(detection = ifelse(station %in% site.extra$Station & detection !=0, 3, detection))

write.csv(df.pred, "SurveyDataWithOffsets.csv", row.names=FALSE)

#16. Check # of recordings per site----
df.size <- data.frame(table(df.pred$station, df.pred$year, df.pred$survey)) %>% 
  rename(station = Var1, year = Var2, survey = Var3) %>% 
  dplyr::filter(survey=="ARU", !Freq %in% c(0, 30))
df.size
