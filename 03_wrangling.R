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

#1. Read in data & wrangle----
files <- fs::dir_ls("/Users/ellyknight/Documents/Employment/Contracts/WatertonCONIMonitoring/Processing/Results/Validated/", glob = "*validated*.txt")

df.raw <- readr::read_tsv(files, 
                      col_names = c("Filename", "TimeOffset", "Duration", "Level", "Quality", "Score", "Recognizer", "Comments"),
                      col_types = "ccddddcc")

#2. Identify non-detections from raw validated data----
df.0 <- df.raw %>% 
  separate(Filename, 
           into = c("Drive", "Park","Batch", "Station", "Date", "Time", "Filetype"),
           sep = "[^A-Za-z0-9 \\- :]",
           remove = FALSE) %>% # NOT letters, numbers, dash, or colon
  select(-Drive, -Park, -Batch, -Filetype, -Recognizer, -Duration, -Level, -Quality, -Score) %>% 
  mutate(Date = ymd(Date),
         DateTime = as_datetime(paste0(Date, Time)), .after = Time) %>% 
  rename(detection=Comments) %>% 
  dplyr::filter(detection==0)

#3. Read in detections & tidy----
df.1 <- readr::read_tsv("/Users/ellyknight/Documents/Employment/Contracts/WatertonCONIMonitoring/Processing/Results/Validated/TruePos_only_ECK.txt",
                        col_names = c("Filename", "TimeOffset", "Duration", "Level", "Quality", "Score", "Recognizer", "Comments"),
                        col_types = "ccddddcc") %>% 
  separate(Filename, 
           into = c("Blank", "Volumes", "Drive", "Park","Batch", "Station", "Date", "Time", "Filetype"),
           sep = "[^A-Za-z0-9 \\- :]",
           remove = FALSE) %>% # NOT letters, numbers, dash, or colon
  select(-Blank, -Volumes, -Drive, -Park, -Batch, -Filetype, -Recognizer, -Duration, -Level, -Quality, -Score) %>% 
  mutate(Date = ymd(Date),
         DateTime = as_datetime(paste0(Date, Time)), .after = Time) %>% 
  rename(detection=Comments)

#4. Put validated data together----
df <- rbind(df.0, df.1)

write.csv(df, "Data/ValidatedRecognizerResults.csv", row.names=FALSE)
