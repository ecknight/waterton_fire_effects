library(tidyverse)
library(sf)
library(raster)

#1. Study area----
dat <- read.csv("SurveyDataWithOffsets&Covariates.csv") %>% 
  mutate(Kenow = ifelse(FireHistory==2017, "impact", "control"),
         Kenow = ifelse(is.na(Kenow), "control", Kenow),
         date = date(DateTime)) %>% 
  dplyr::select(survey, station, year, X, Y) %>% 
  unique()
