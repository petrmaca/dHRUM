library(dHRUM)
library(data.table)
library(fst)

SRAdta <- as.data.table(read.fst("/home/hubert/MyData/POH/data/OUT/chp4Plus/chp4Plus_SRA_SC2.fst"))
TEMPMindta <- as.data.table(read.fst("/home/hubert/MyData/POH/data/OUT/chp4Plus/chp4Plus_TMIN_SC2.fst"))
TEMPMaxdta <- as.data.table(read.fst("/home/hubert/MyData/POH/data/OUT/chp4Plus/chp4Plus_TMAX_SC2.fst"))
TEMPavgdta <- merge(TEMPMindta, TEMPMaxdta, by = c("DTM", "chp_4"))
TEMPavgdta[,Tavg := (Tmin + Tmax)/2]
Qmd <-  as.data.table(read.fst("/home/hubert/MyData/POH/data/rmvody/mQ_4chpPlus.fst"))

nHrus <- 1
ntreads <- 1
#Areas <- runif(nHrus,min = 1,max  = 10) #[m2]
Areas <- runif(nHrus,min = 38780000,max  = 38780050)
IdsHrus <- paste0("ID",seq(1:length(Areas)))

