library(tidyverse)
library(suncalc)
library(lubridate)

#Get server files
files <- list.files("/Volumes/BUpublic/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE)

files.df <- data.frame(path=files) %>% 
  separate(path, into=c("f1", "f2", "f3", "project", "year", "site", "AB", "f4", "recording"), sep="/", remove=FALSE) %>% 
  mutate(recording=ifelse(AB %in% c("A", "B"), recording, f4),
         recording=ifelse(str_sub(AB, -4, -1) %in% c(".wav", ".wac"), AB, recording),
         round=NA,
         station=NA) %>% 
  dplyr::filter(!is.na(recording),
                str_sub(recording, -4, -1) %in% c(".wav", ".wac")) %>% 
  dplyr::select(path, project, year, round, station, site, recording)

#NOTE WLNP-017-03 2021 files are misnamed

#write.csv(files.df, "ServerRecordingList_2021.csv", row.names = FALSE)

files.server <- rbind(read.csv("ServerRecordingList.csv"), read.csv("ServerRecordingList_2021.csv")) %>% 
  mutate(location="Server")

#Get hard drive files
files.2020 <- list.files("/Volumes/Contracts/4 - 2020 WLNP ARU Pilot/2 - Raw Acoustic Data/", recursive=TRUE, full.names = TRUE, include.dirs=TRUE)

files.df <- data.frame(path=files.2020) %>% 
  separate(path, into=c("f1", "f2", "f3", "f4", "f5", "f6", "station", "site", "f7", "recording"), sep="/", remove=FALSE) %>% 
  mutate(year="2020", round=NA, project="WLNP") %>% 
  dplyr::filter(!is.na(recording),
                str_sub(recording, -4, -1)==".wav") %>% 
  dplyr::select(path, project, year, round, station, site, recording)

#write.csv(files.df, "HardDriveRecordingList_2020.csv", row.names = FALSE)

files.2021 <- list.files("/Volumes/LaCie/5 - 2021 WLNP ARU/2 - Raw Acoustic Data", recursive=TRUE, full.names = TRUE, include.dirs=TRUE)

files.df <- data.frame(path=files.2021) %>% 
  separate(path, into=c("f1", "f2", "f3", "f4", "f5", "station", "site", "f6", "recording"), sep="/", remove=FALSE) %>% 
  mutate(year="2021", round=NA, project="WLNP") %>% 
  dplyr::filter(!is.na(recording),
                str_sub(recording, -4, -1)==".wav") %>% 
  dplyr::select(path, project, year, round, station, site, recording)

#write.csv(files.df, "HardDriveRecordingList_2021.csv", row.names = FALSE)

files.hd <- rbind(read.csv("HardDriveRecordingList_2020.csv"), read.csv("HardDriveRecordingList_2021.csv")) %>% 
  mutate(location="HD")

files.all <- rbind(files.server, files.hd)

#Get time relative to sunset
#Use mean X Y for existing locations
locs <- read.csv("ExistingSamplingLocations.csv") %>% 
  dplyr::select(X, Y) %>% 
  summarize(X=mean(X),
            Y=mean(Y))

files.locs <- files.all %>% 
  mutate(lon=locs$X,
         lat=locs$Y) %>% 
  separate(recording, into=c("site2", "datename", "timename"), sep="_", remove=FALSE) %>% 
  mutate(date=ymd(paste(str_sub(datename, 1, 4),
                        str_sub(datename, 5, 6),
                        str_sub(datename, 7, 8),
                        sep="-")),
         hour=str_sub(timename, 1, 2),
         min=str_sub(timename, 3, 4),
         yday=yday(date))

files.sun <- getSunlightTimes(data=files.locs, tz="Canada/Mountain") %>% 
  cbind(files.locs %>% 
          dplyr::select(-lat, -lon, -date)) %>% 
  mutate(datetime = ymd_hm(paste0(as.character(date), " ", hour, ":", min), tz="Canada/Mountain"),
         tsss = as.numeric(difftime(datetime, sunset, units="hours")),
         tssr = as.numeric(difftime(datetime, sunrise, units="hours")),
         tsss = ifelse(tsss < -12, tsss+24, tsss),
         tssr = ifelse(tssr > 12, tssr-24, tssr)) 

#Inspect temporal distribution
ggplot(files.sun) +
  geom_histogram((aes(x=tsss))) +
  facet_wrap(location~year)

ggplot(files.sun) +
  geom_histogram((aes(x=yday))) +
  facet_wrap(location~year)

#Filter by suntimes & date
files.p <- files.sun %>% 
  dplyr::filter(tsss >= -2, tsss <= 6,
                yday >= 152, yday <= 213)

table(files.p$site, files.p$year)

#Filter out sites without enough data
N <- 30

files.n <- data.frame(table(files.p$site, files.p$year)) %>% 
  rename(site=Var1, year=Var2) %>% 
  dplyr::filter(Freq > N) %>% 
  mutate(year=as.numeric(as.character(year)))

files.total <- length(unique(files.sun$site))

#Random sample
set.seed(1234)
files.select <- files.p %>% 
  right_join(files.n) %>% 
  group_by(year, site) %>% 
  sample_n(30) %>% 
  ungroup()

table(files.select$site, files.select$year)

ggplot(files.select) +
  geom_histogram((aes(x=tsss))) +
  facet_wrap(location~year)

ggplot(files.select) +
  geom_histogram((aes(x=yday))) +
  facet_wrap(location~year)

#Set up folders for processing
files.copy <- files.select %>% 
  mutate(file=ceiling(row_number()/300))

for(i in 1:max(files.copy$file)){
  if(!dir.exists(file.path(paste0("/Volumes/Seagate4TB/WLNP/", i)))){
    dir.create(file.path(paste0("/Volumes/Seagate4TB/WLNP/", i)))
  }
}

#Copy
for(i in 4620:nrow(files.copy)){
  from.i <- files.copy$path[i]
  to.i <- paste0("/Volumes/Seagate4TB/WLNP/", files.copy$file[i], "/", files.copy$recording[i])
  file.copy(from.i, to.i)
  print(paste0("Completed copying recording ", i, " of ", nrow(files.copy), " recordings"))
}

#Check 
files.0 <- file.info(list.files("/Volumes/Seagate4TB/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE)) %>%
  dplyr::filter(size==0)
files.0$path <- row.names(files.0)
files.0.sep <- files.0 %>% 
  separate(path, into=c("f1", "f2", "f3", "f4", "folder", "file"), sep="/", remove=FALSE) %>% 
  separate(file, into=c("site", "datename", "timename"), sep="_", remove=FALSE) %>% 
  mutate(year=str_sub(datename, 1, 4))

#Select new files
files.new <- data.frame(table(files.0.sep$site, files.0.sep$year, files.0.sep$folder)) %>% 
  rename(site = Var1, year=Var2, folder=Var3) %>% 
  dplyr::filter(Freq > 0)

files.1 <- file.info(list.files("/Volumes/Seagate4TB/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE)) %>%
  dplyr::filter(size!=0)
files.1$path <- row.names(files.1)
files.1.sep <- files.1 %>% 
  separate(path, into=c("f1", "f2", "f3", "f4", "folder", "file"), sep="/", remove=FALSE) %>% 
  separate(file, into=c("site", "datename", "timename"), sep="_", remove=FALSE) %>% 
  dplyr::filter(!is.na(timename))
files.p.available <- files.p %>% 
  dplyr::filter(!recording %in% files.1.sep$file)

set.seed(1234)
files.select2 <- data.frame()
for(i in 1:nrow(files.new)){
  
  files.p.i <- files.p.available %>% 
    dplyr::filter(site==files.new$site[i],
                  year==files.new$year[i])
  
  files.select2.i <- files.p.i %>% 
    sample_n(files.new$Freq[i]) %>% 
    dplyr::select(site, year, path, recording) %>% 
    mutate(folder=files.new$folder[i])
  
  files.select2 <- rbind(files.select2, files.select2.i)
  
}

#Copy
for(i in 1:nrow(files.select2)){
  from.i <- files.select2$path[i]
  to.i <- paste0("/Volumes/Seagate4TB/WLNP/", files.select2$folder[i], "/", files.select2$recording[i])
  file.copy(from.i, to.i)
  print(paste0("Completed copying recording ", i, " of ", nrow(files.select2), " recordings"))
}

#Remove empty files
file.remove(files.0$path)

#Check that all sites have 30 recordings
files.use <- file.info(list.files("/Volumes/Seagate4TB/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE))
files.use$path <- row.names(files.use)
files.use.sep <- files.use %>% 
  separate(path, into=c("f1", "f2", "f3", "f4", "folder", "recording"), sep="/", remove=FALSE) %>% 
  separate(recording, into=c("site", "datename", "timename"), sep="_", remove=FALSE) %>% 
  mutate(year=as.numeric(str_sub(datename, 1, 4))) %>% 
  dplyr::filter(!is.na(timename)) %>% 
  dplyr::select(site, year, recording)

files.count <- data.frame(table(files.use.sep$site, files.use.sep$year)) %>% 
  rename(site=Var1, year=Var2) %>% 
  dplyr::filter(Freq < 30, Freq > 0)

#Select new files
files.1 <- file.info(list.files("/Volumes/Seagate4TB/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE)) %>%
  dplyr::filter(size!=0)
files.1$path <- row.names(files.1)
files.1.sep <- files.1 %>% 
  separate(path, into=c("f1", "f2", "f3", "f4", "folder", "file"), sep="/", remove=FALSE) %>% 
  separate(file, into=c("site", "datename", "timename"), sep="_", remove=FALSE) %>% 
  dplyr::filter(!is.na(timename))
files.p.available <- files.p %>% 
  dplyr::filter(!recording %in% files.1.sep$file)

set.seed(1234)
files.select3 <- data.frame()
for(i in 3:nrow(files.count)){
  
  files.p.i <- files.p.available %>% 
    dplyr::filter(site==files.count$site[i],
                  year==files.count$year[i])
  
  files.select3.i <- files.p.i %>% 
    sample_n(30 - files.count$Freq[i]) %>% 
    dplyr::select(site, year, path, recording) %>% 
    mutate(folder=files.count$folder[i])
  
  files.select3 <- rbind(files.select3, files.select3.i)
  
}

#Save out list of files
files.all <- file.info(list.files("/Volumes/Seagate4TB/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE))
files.all$path <- row.names(files.all)

write.csv(files.all, "FinalFileList.csv", row.names = FALSE)

#THESE FOUR SITES HAVE NO MORE RECORDINGS. CONSIDER REMOVING

#Check filenames of 2021 data
files.2021 <- read.csv("ServerRecordingList_2021.csv") %>% 
  separate(recording, into=c("site2", "datename", "timename"), sep="_", remove=FALSE) %>% 
  dplyr::filter(site!=site2)

files.use <- file.info(list.files("/Volumes/Seagate4TB/WLNP", recursive=TRUE, full.names = TRUE, include.dirs=TRUE))
files.use$path <- row.names(files.use)
files.use.sep <- files.use %>% 
  separate(path, into=c("f1", "f2", "f3", "f4", "folder", "recording"), sep="/", remove=FALSE) %>% 
  separate(recording, into=c("site", "datename", "timename"), sep="_", remove=FALSE) %>% 
  mutate(year=as.numeric(str_sub(datename, 1, 4))) %>% 
  dplyr::filter(!is.na(timename)) %>% 
  dplyr::select(site, year, recording)

files.names <- files.2021 %>% 
  inner_join(files.use.sep) 