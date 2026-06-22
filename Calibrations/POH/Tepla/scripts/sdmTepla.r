library(dHRUM)
library(data.table)
library(fst)
library(lubridate)

Qmd <-  as.data.table(read.fst("/home/hubert/MyData/POH/data/rmvody/mQ_4chpPlus.fst"))
dta <- as.data.table(readRDS("/home/hubert/MyData/POH/data/Tepla/IN/INP_Tepla.rds"))
ddd <- as.data.table(readRDS("/home/hubert/MyData/POH/data/Tepla/IN/INF_Tepla.rds"))
ddd[,ID := chp_14_s]

teplaCHP <- as.data.table(read.table("/home/hubert/MyData/POH/data/OUT/catchments/Tepla_chp14s.txt",header = 1))

chps_Brezova <- teplaCHP[Povodi %in% "Brezova",]
chps_Brezova[, ID := chp14s]

resBrezova <- merge(dta,chps_Brezova, by = "ID", all.y = TRUE )

resBrezovaII <- merge(ddd,resBrezova, by = "ID", all.y = TRUE )

length(resBrezovaII[, unique(ID)])
nrow(chps_Brezova)

Days <- c(30,60,90,120,150,180,210,240,270,300,330,355,364)

