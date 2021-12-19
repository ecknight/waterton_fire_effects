library(tidyverse)
library(tidylog)

#Note: This file is from the bioacoustics score manuscript processing

#1. Read in & wrangle data for a score of 60----
dist <- read.csv("Data/FinalClipListAndRecognizerResults.csv")

dist.1 <- dist %>% 
  dplyr::filter(use!="N",
                score >= 60, 
                recognizer == "near training data") %>% 
  dplyr::select(distance, file, clip, score) %>% 
  mutate(detection = 1)

dist.all <- dist %>% 
  dplyr::filter(use!="N") %>% 
  dplyr::select(distance, file, clip) %>% 
  unique() %>% 
  left_join(dist.1) %>% 
  mutate(detection = ifelse(is.na(detection), 0, detection))
table(dist.all$detection)

#2. Fit logist regression----
md1 <- glm(detection ~ I(-distance^2) - 1, dist.all, family=binomial("log"))
md2 <- glm(detection ~ I(-distance^2) - 1, dist.all, family=binomial("cloglog"))
AIC(md1, md2)

#3. Calculate EDR----
#' We pick `md1` over `md2` based on AIC
EDR <- data.frame(edr = unname(sqrt(1/coef(md1)))) %>% 
  mutate(eda = edr^2*pi*0.0001)
EDR
#       edr      eda
#1 247.7956 19.29021

write.csv(EDR, "EDR.csv", row.names = FALSE)

#4. Calculate recall----
sum(dist.all$detection)/nrow(dist.all)
#0.5626959