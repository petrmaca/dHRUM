library(dHRUM)
library(data.table)
library(fst)
library(lubridate)

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

Days <- c(30,60,90,120,150,180,210,240,270,300,330,355,364)

basinCHP <- Qmd[1,chp_4]

sra <- SRAdta[chp_4 %in% basinCHP,SRA]
temp <- TEMPavgdta[chp_4 %in% basinCHP,Tavg]

LAdt = as.data.table(get_Lai_Data())
ctgrs = unique(LAdt$Category)
ctgrs
lai = LAdt[Category %in% ctgrs[17],mean_LAI]
ll = nrow(SRAdta[chp_4 %in% basinCHP,])
ntms <- round(ll/364.25) + 1

laiALL <- rep(x = lai, times = 61)

DTall <- data.table(sra = sra, temp = temp, lai = laiALL[1:ll])


nHrus <- 1
ntreads <- 1
#Areas <- runif(nHrus,min = 1,max  = 10) #[m2]
Areas <- runif(nHrus,min = 38780000,max  = 38780050)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus,ntreads)

setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep("DDF",times = length(Areas)),hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("Pitman",times= length(Areas)))
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)

pars <-getCurdHRUpars(dHRUM_ptr = dhrus, singleHruId = 0)

