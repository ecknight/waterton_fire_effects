library(tidyverse)
library(lubridate)

#1. Read in data and wrangle----
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

#2. Identify recordings with wingbooms----
rec.b <- val %>% 
  dplyr::filter(Comments=="B") %>% 
  dplyr::select(Site, Filename) %>% 
  unique() %>% 
  group_by(Site) %>% 
  sample_n(1) %>% 
  ungroup()

#3. Identify recordings with calls----
rec.c <- val %>% 
  dplyr::select(Site, Filename) %>% 
  unique() %>% 
  group_by(Site) %>% 
  sample_n(1) %>% 
  ungroup()

#4. Bind to validation results----
val.b <- val %>% 
  dplyr::filter(Filename %in% rec.b$Filename) %>% 
  mutate(Individual="") %>% 
  dplyr::select("Filename", "TimeOffset", "Duration", "Level", "Quality", "Score", "Comments", "Individual")

val.c <- val %>% 
  dplyr::filter(Filename %in% rec.c$Filename) %>% 
  mutate(Individual="") %>% 
  dplyr::select("Filename", "TimeOffset", "Duration", "Level", "Quality", "Score", "Comments", "Individual")

#4. Write out----
write.table(val.b, "IndividualID_boom.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(val.c, "IndividualID_call.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)

