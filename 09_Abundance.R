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

#4. Sandwich estimator for variance----
Xp <- X1[!duplicated(abu.boom$station),]
rownames(Xp) <- abu$sa[!duplicated(abu$sa)]

#Boom
CF.boom <- NULL
for (i in 1:300) {
  abux <- NULL
  abuj <- abu.boom[sample.int(nrow(abu.boom), replace=TRUE),]
  abux <- rbind(abux, abuj[!duplicated(abuj$station),])
  
  Y1i <- abux$abundance
  X1i <- model.matrix(~ 1, abux)
  
  cl1pi <- try(zi.fit(Y1i, X1i, distr="pois", type="CL", hessian=TRUE)$CL)
  if (!inherits(cl1pi, "try-error"))
    est <- data.frame(lambda = exp(cl1pi$coef),
                      i = i)
  CF.boom <- rbind(CF.boom, est)
}

#call
CF.call <- NULL
for (i in 1:300) {
  abux <- NULL
  abuj <- abu.call[sample.int(nrow(abu.call), replace=TRUE),]
  abux <- rbind(abux, abuj[!duplicated(abuj$station),])
  
  Y1i <- abux$abundance
  X1i <- model.matrix(~ 1, abux)
  
  cl1pi <- try(zi.fit(Y1i, X1i, distr="pois", type="CL", hessian=TRUE)$CL)
  if (!inherits(cl1pi, "try-error"))
    est <- data.frame(lambda = exp(cl1pi$coef),
                      i = i)
  CF.call <- rbind(CF.call, est)
}

#extraterritorial
CF.extra <- NULL
for (i in 1:300) {
  abux <- NULL
  abuj <- abu.extra[sample.int(nrow(abu.extra), replace=TRUE),]
  abux <- rbind(abux, abuj[!duplicated(abuj$station),])
  
  Y1i <- abux$abundance
  X1i <- model.matrix(~ 1, abux)
  
  cl1pi <- try(zi.fit(Y1i, X1i, distr="pois", type="CL", hessian=TRUE)$CL)
  if (!inherits(cl1pi, "try-error"))
    est <- data.frame(lambda = exp(cl1pi$coef),
                      i = i)
  CF.extra <- rbind(CF.extra, est)
}

#5. Estimates----

#Boom
# mean of the poisson (including 0 counts and offsets too)
lambda.boom <- exp(cl0p.boom$coef)
lambda.boom
# mean of the sandwich estimator
lambda.boom <- CF.boom %>% 
  summarize(quantlow = quantile(lambda, prob=0.025),
            quanthigh = quantile(lambda, 0.975),
            lambda=mean(lambda)) %>% 
  mutate(detection = "boom")
lambda.boom
# P(N=0) based on Poisson, which is the suitable but unoccupied probability
exp(-lambda.boom$lambda)
## suitable and occupied
1-exp(-lambda.boom$lambda)

#Call
# mean of the poisson (including 0 counts and offsets too)
lambda.call <- exp(cl0p.call$coef)
lambda.call
## mean of sandwich estimate
lambda.call <- CF.call %>% 
  summarize(quantlow = quantile(lambda, prob=0.025),
            quanthigh = quantile(lambda, 0.975),
            lambda=mean(lambda)) %>% 
  mutate(detection = "call")
lambda.call
# P(N=0) based on Poisson, which is the suitable but unoccupied probability
exp(-lambda.call$lambda)
# suitable and occupied
1-exp(-lambda.call$lambda)

#Extraterritorial
## mean of the poisson (including 0 counts and offsets too)
lambda.extra <- exp(cl0p.extra$coef)
lambda.extra
## mean of sandwich estimate
lambda.extra <- CF.extra %>% 
  summarize(quantlow = quantile(lambda, prob=0.025),
            quanthigh = quantile(lambda, 0.975),
            lambda=mean(lambda)) %>% 
  mutate(detection = "extra")
lambda.extra
## P(N=0) based on Poisson, which is the suitable but unoccupied probability
exp(-lambda.extra$lambda)
## suitable and occupied
1-exp(-lambda.extra$lambda)

#6. Save out----
lambda <- rbind(lambda.boom, lambda.call, lambda.extra)
lambda

write.csv(lambda, "LambdaEstimates.csv", row.names = FALSE)
