library(tidyverse)

locs.2020 <- read.csv("CommonNightawkSamplingLocationsExisting.csv") %>% 
  rename(Latitude=Y, Longitude=X)
locs.2021 <- read.csv("CommonNightawkSamplingLocations2021.csv") %>% 
  dplyr::select(Site, Latitude, Longitude)
locs.aru <- rbind(locs.2020, locs.2021)

locs.all <- read.csv("ExistingSamplingLocations.csv") %>% 
  rename(PresenceBefore = Presence,
         Longitude = X,
         Latitude = Y) 

recs <- read.csv("FinalFileList.csv")

val <- readr::read_tsv("TruePos_only_ECK.txt", 
                       col_names = c("Filename", "TimeOffset", "Duration", "Level", "Quality", "Score", "Recognizer", "Comments"),
                       col_types = "ccddddcc") %>% 
  separate(Filename, 
           into = c("Drive1", "Drive2", "Drive3", "Park","Batch", "Site", "Date", "Time", "Filetype"),
           sep = "[^A-Za-z0-9 \\- :]",
           remove = FALSE) %>% # NOT letters, numbers, dash, or colon
  select(-Drive1, -Drive2, -Drive3, -Park, -Batch, -Filetype) %>% 
  mutate(Date = ymd(Date),
         DateTime = as_datetime(paste0(Date, Time)), .after = Time,
         Year=year(DateTime)) %>% 
  dplyr::filter(Comments!="0")

table(val$Comments)

val.sum <- val %>% 
  group_by(Site, Year, Comments, Filename) %>% 
  summarize(n=n()) %>% 
  mutate(pres = ifelse(n>0, 1, 0)) %>% 
  group_by(Site, Year, Comments) %>% 
  summarize(recs = sum(pres)) %>% 
  ungroup()
left_join(locs.aru %>% 
            mutate(Site = gsub(pattern="CONI-GAP", replacement="GAP", x=Site),
                   Site = gsub(pattern="Cardston Gate 0", replacement = "CAR-", x=Site),
                   Site = gsub(pattern="Cardson Entrance 0", replacement = "ENT-", x=Site),
                   Site = gsub(pattern="Eskerine 0", replacement = "ESK-", x=Site),
                   Site = gsub(pattern="HORSESHOE-", replacement = "HORSE-", x=Site),
                   Site = gsub(pattern="LAKEVIEW-", replacement = "LAKE-", x=Site),
                   Site = gsub(pattern="OILBASIN", replacement = "OIL", x=Site),
                   Site = gsub(pattern="Red Rock Parkway ", replacement = "RRP-", x=Site),
                   Site = gsub(pattern="WLNP-006-001", replacement = "WLNP-6-1", x=Site),
                   Site = gsub(pattern="WLNP-007-001", replacement = "WLNP-7-1", x=Site)))

table(val.sum$Comments)

val.na <- val.sum %>% 
  dplyr::filter(is.na(Latitude))


