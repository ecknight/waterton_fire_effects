#title: Availability for detection offsets for CONI density modelling
#author: Elly C. Knight, adapted from Peter Solymos
#created: October 31, 2021

library(mefa4)
library(maptools)
library(survival)
library(tidyverse)
library(readr)
library(fs)
library(lubridate)

options(scipen=99999)

#This step uses just the new recognizer data to model probability of availability for date and time and then uses it to predict offsets for human and recognizer data that will be used in the models

#1. Read in recognizer data----
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
                !Filename %in% df.raw$Station) %>% 
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

#3. Add new fields-----

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
  dplyr::select(df, loglik, AICc, delta, weight)
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
aic.call

#Seems like cos + sin fits best for both

#10. Summary of predictions----
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

#11. Visualize predictions----

#Boom
vjd.boom <- seq(min(df.boom$jday), max(df.boom$jday), len=51*10)
vtd.boom <- seq(0, 23/24, len=24*10)
pr.boom <- expand.grid(jday=vjd.boom, tod=vtd.boom)
pr.boom$jday2 <- pr.boom$jday^2
pr.boom$cos <- cos(pr.boom$tod * 2 * pi)
pr.boom$sin <- sin(pr.boom$tod * 2 * pi)
#' make predictions for the grid
fit.boom <- 1/predict(mb.boom, newdata=pr.boom)
summary(fit.boom)
summary(1-exp(-fit.boom*10))
#' we plot contours for P(detection in 10 min)
z.boom <- matrix(1-exp(-fit.boom*10), length(vjd.boom), length(vtd.boom))
plot(tod*24 ~ jitter(doy), df.boom, pch=19, cex=0.6, col="#80808020",
     main="Probability of availability (10 min)",
     xlab="Ordinal day", ylab="Hour")
points(tod*24 ~ doy, df.boom[df.boom$detection > 0,], pch=19, cex=0.6,
       col="#FF0040")
l <- c(0.001, 0.01, 0.1, 0.2, 0.4)
for (i in seq_along(l)) {
  col <- colorRampPalette(c("blue", "red"))( length(l) )[i]
  contour(vjd.boom*365, vtd.boom*24, z.boom, add=TRUE, col=col, levels=l[i],
          labcex=0.8, lwd=1)
}

#Call
vjd.call <- seq(min(df.call$jday), max(df.call$jday), len=51*10)
vtd.call <- seq(0, 23/24, len=24*10)
pr.call <- expand.grid(jday=vjd.call, tod=vtd.call)
pr.call$jday2 <- pr.call$jday^2
pr.call$cos <- cos(pr.call$tod * 2 * pi)
pr.call$sin <- sin(pr.call$tod * 2 * pi)
#' make predictions for the grid
fit.call <- 1/predict(mb.call, newdata=pr.call)
summary(fit.call)
summary(1-exp(-fit.call*10))
#' we plot contours for P(detection in 10 min)
z.call <- matrix(1-exp(-fit.call*10), length(vjd.call), length(vtd.call))
plot(tod*24 ~ jitter(doy), df.call, pch=19, cex=0.6, col="#80808020",
     main="Probability of availability (10 min)",
     xlab="Ordinal day", ylab="Hour")
points(tod*24 ~ doy, df.call[df.call$detection > 0,], pch=19, cex=0.6,
       col="#FF0040")
l <- c(0.001, 0.01, 0.1, 0.2, 0.4)
for (i in seq_along(l)) {
  col <- colorRampPalette(c("blue", "red"))( length(l) )[i]
  contour(vjd.call*365, vtd.call*24, z.call, add=TRUE, col=col, levels=l[i],
          labcex=0.8, lwd=1)
}

#12. Calculate offsets for recognizer data----
df.aru <- df %>% 
  mutate(detection = ifelse(detection=="B", 2, as.numeric(detection))) %>% 
  group_by(Filename, Station, DateTime) %>% 
  summarize(detection=max(detection)) %>% 
  ungroup() %>% 
  mutate(doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         tod = start/24,
         sin = sin(tod*2*pi),
         cos = cos(tod*2*pi),
         survey="ARU") %>% 
  rename(station=Station)
df.aru$p.boom <- 1-exp(-(1/predict(mb.boom, newdata = df.aru))*10)
df.aru$p.call <- 1-exp(-(1/predict(mb.call, newdata = df.aru))*10)

#13. Calculate offsets for human data----
df.hum <- read.csv("Data/HumanData.csv") %>% 
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
         doy = yday(DateTime),
         jday = doy/365,
         jday2 = jday^2,
         start = hour(DateTime) + minute(DateTime)/60,
         tod = start/24,
         sin = sin(tod*2*pi),
         cos = cos(tod*2*pi),
         survey="human") %>% 
  rename(station=Site, detection=Presence)
df.hum$p.boom <- 1-exp(-(1/predict(mb.boom, newdata = df.hum))*6)
df.hum$p.call <- 1-exp(-(1/predict(mb.call, newdata = df.hum))*6)


#14. Correct for differences in recorder & recognizer EDR----
#From Yip et al. Fig 2A SM2 @ 4000 Hz
SM2 <- 0.7 

#Recall at score 60 from known distance clips
rec <-  0.5626959

df.hum$p.boom <- df.hum$p.boom*(1/SM2)*(1/rec)
df.hum$p.call <- df.hum$p.call*(1/SM2)*(1/rec)

#15. Put together----
df.pred <- rbind(df.aru %>% 
                   dplyr::select(survey, station, DateTime, detection, p.boom, p.call),
                 df.hum  %>% 
                   dplyr::select(survey, station, DateTime, detection, p.boom, p.call)) %>% 
  mutate(year=year(DateTime),
         fire=ifelse(year <=2017, "before", "after"))

write.csv(df.pred, "SurveyDataWithOffsets.csv", row.names=FALSE)
