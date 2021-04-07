library(tidyverse)

files <- list.files("/Volumes/BUpublic/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE)

files.df <- data.frame(path=files) %>% 
  separate(path, into=c("f1", "f2", "f3", "project", "year", "round", "station", "site", "recording"), sep="/", remove=FALSE) %>% 
  dplyr::filter(!is.na(recording),
                str_sub(recording, -4, -1)==".wav") %>% 
  dplyr::select(-f1, -f2, -f3)

write.csv(files.df, "ServerRecordingList.csv", row.names = FALSE)
