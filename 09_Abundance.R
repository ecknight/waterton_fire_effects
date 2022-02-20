#title: ZIP abundance model fitting for CONI density modelling
#author: Elly C. Knight, adapted from Peter Solymos
#created: November 8, 2021

library(tidyverse)
library(lubridate)

options(scipen=99999)

source("src/zi.fit.R")

#1. Wrangling----
dat.aru.b <- read.table("Data/IndividualID_boom_validated.txt")
dat.aru.c <- read.table("Data/IndividualID_call_validated.txt")
colnames(dat.aru.b) <- c("filename", "time", "duration", "level", "quality", "score", "detection", "individual")
colnames(dat.aru.c) <- c("filename", "time", "duration", "level", "quality", "score", "detection", "individual")

abu.aru.b <- dat.aru.b %>% 
  dplyr::filter(individual!="N",
                detection=="B") %>% 
  mutate(individualUse = as.numeric(str_sub(individual, 1, 1))) %>% 
  group_by(filename) %>% 
  summarize(abundance = max(individualUse)) %>% 
  ungroup() %>% 
  mutate(survey="ARU") %>% 
  separate(filename, into=c("blank", "volumes", "hd", "project", "folder", "file"), sep="/", remove=FALSE) %>% 
  separate(file, into=c("station", "datename", "timename"), sep="_", remove=FALSE) %>% 
  mutate(date = ymd(datename),
         year = year(date)) %>% 
  dplyr::select(survey, station, year, date, abundance)

abu.aru.c <- dat.aru.c %>% 
  dplyr::filter(individual!="N",
                detection=="1") %>% 
  mutate(individualUse = as.numeric(str_sub(individual, 1, 1))) %>% 
  group_by(filename) %>% 
  summarize(abundance = max(individualUse)) %>% 
  ungroup() %>% 
  mutate(survey="ARU") %>% 
  separate(filename, into=c("blank", "volumes", "hd", "project", "folder", "file"), sep="/", remove=FALSE) %>% 
  separate(file, into=c("station", "datename", "timename"), sep="_", remove=FALSE) %>% 
  mutate(date = ymd(datename),
         year = year(date)) %>% 
  dplyr::select(survey, station, year, date, abundance)

abu.aru.x <- abu.aru.c %>% 
  anti_join(abu.aru.b %>% 
              dplyr::select(survey, station, year))

dat.hum <- read.csv("Data/HumanData.csv")

abu.hum <- dat.hum %>% 
  rename(abundance = X..Individuals, survey=SurveyType, date=Date, year=Year, station=Site) %>% 
  dplyr::filter(abundance > 0) %>% 
  mutate(detection = "call") %>% 
  dplyr::select(survey, station, year, date, detection, abundance)

abu.boom <- abu.aru.b %>% 
  mutate(fire = ifelse(year <= 2017, "before", "after")) %>% 
  mutate(firesurvey = paste0(fire,"-", survey))

abu.call <- abu.aru.c %>% 
  mutate(fire = ifelse(year <= 2017, "before", "after")) %>% 
  mutate(firesurvey = paste0(fire,"-", survey))

abu.extra <- abu.aru.c %>% 
  anti_join(abu.aru.b %>% 
              dplyr::select(survey, station, year)) %>% 
  mutate(fire = ifelse(year <= 2017, "before", "after")) %>% 
  mutate(firesurvey = paste0(fire,"-", survey))
  

#2. Visualize----
ggplot(abu.boom) +
  geom_histogram(aes(x=abundance)) +
  facet_grid(fire ~ survey)

ggplot(abu.call) +
  geom_histogram(aes(x=abundance)) +
  facet_grid(fire ~ survey)

ggplot(abu.extra) +
  geom_histogram(aes(x=abundance)) +
  facet_grid(fire ~ survey)

#3. Model----

#Boom
Y1 <- abu.boom$abundance
X0 <- data.matrix(rep(1, length(Y1)))
Z0 <- X0

cl0p <- zi.fit(Y1, X0, Z0, distr="pois", type="CL", hessian=TRUE)$CL
cl0nb <- zi.fit(Y1, X0, Z0, distr="negbin", type="CL", hessian=TRUE)$CL

ic <- AIC(cl0p, cl0nb)
ic$loglik <- c(logLik(cl0p), logLik(cl0nb))
ic$BIC <- AIC(cl0p, cl0nb, k=log(length(Y1)))$AIC
ic$AICc <- ic$AIC + (2*ic$df^2+2*ic$df) / (length(Y1)-ic$df-1)
ic <- ic[order(ic$AICc),] %>% 
  mutate(delta = AICc - head(AICc, 1),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>% 
  dplyr::select(df, loglik, AICc, delta, weight)
ic

cl0p.boom <- cl0p


#Call
Y1 <- abu.call$abundance
X0 <- data.matrix(rep(1, length(Y1)))
Z0 <- X0

cl0p <- zi.fit(Y1, X0, Z0, distr="pois", type="CL", hessian=TRUE)$CL
cl0nb <- zi.fit(Y1, X0, Z0, distr="negbin", type="CL", hessian=TRUE)$CL

ic <- AIC(cl0p, cl0nb)
ic$loglik <- c(logLik(cl0p), logLik(cl0nb))
ic$BIC <- AIC(cl0p,cl0nb, k=log(length(Y1)))$AIC
ic$AICc <- ic$AIC + (2*ic$df^2+2*ic$df) / (length(Y1)-ic$df-1)
ic <- ic[order(ic$AICc),] %>% 
  mutate(delta = AICc - head(AICc, 1),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>% 
  dplyr::select(df, loglik, AICc, delta, weight)
ic

cl0p.call <- cl0p

#Extraterritorial
Y1 <- abu.extra$abundance
X0 <- data.matrix(rep(1, length(Y1)))
Z0 <- X0

cl0p <- zi.fit(Y1, X0, Z0, distr="pois", type="CL", hessian=TRUE)$CL
cl0nb <- zi.fit(Y1, X0, Z0, distr="negbin", type="CL", hessian=TRUE)$CL

ic <- AIC(cl0p, cl0nb)
ic$loglik <- c(logLik(cl0p), logLik(cl0nb))
ic$BIC <- AIC(cl0p,cl0nb, k=log(length(Y1)))$AIC
ic$AICc <- ic$AIC + (2*ic$df^2+2*ic$df) / (length(Y1)-ic$df-1)
ic <- ic[order(ic$AICc),] %>% 
  mutate(delta = AICc - head(AICc, 1),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>% 
  dplyr::select(df, loglik, AICc, delta, weight)
ic

cl0p.extra <- cl0p
#Poisson  with no covs for both detection types

#4. Estimate

#Boom
## mean of the poisson (including 0 counts and offsets too)
lambda.boom <- exp(cl0p.boom$coef)
lambda.boom
## P(N=0) based on Poisson, which is the suitable but unoccupied probability
exp(-lambda.boom)
## suitable and occupied
1-exp(-lambda.boom)

#Call
## mean of the poisson (including 0 counts and offsets too)
lambda.call <- exp(cl0p.call$coef)
lambda.call
## P(N=0) based on Poisson, which is the suitable but unoccupied probability
exp(-lambda.call)
## suitable and occupied
1-exp(-lambda.call)

#Extraterritorial
## mean of the poisson (including 0 counts and offsets too)
lambda.extra <- exp(cl0p.extra$coef)
lambda.extra
## P(N=0) based on Poisson, which is the suitable but unoccupied probability
exp(-lambda.extra)
## suitable and occupied
1-exp(-lambda.extra)

#5. Save out----
lambda <- data.frame(detection=c("boom", "call", "extra"),
                     lambda=c(lambda.boom, lambda.call, lambda.extra))
lambda

write.csv(lambda, "LambdaEstimates.csv", row.names = FALSE)
