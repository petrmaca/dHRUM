library(dHRUM)
library(data.table)
library(fst)
library(lubridate)
library(SoilHyP) # SCE-UA

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

nHrus <- 1
ntreads <- 1
IdsHrus <- paste0(basinCHP)
Areas <- resBrezovaII$a_km2[i]
latBas <- resBrezovaII$LAT[i]

i=1
basinCHP <- as.character(resBrezovaII[i,ID])
dhrus <- initdHruModel(nHrus,Areas,IdsHrus,ntreads)

setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep("DDF",times = length(Areas)),hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("Pitman",times= length(Areas)))
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
setPTLInputsToAlldHrus(dhrus, Prec = DTall$sra, Temp = DTall$temp, Lai = DTall$lai, inDate =  as.Date("1961/01/01"))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
