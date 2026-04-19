# library(profvis)
# profvis(
# {
library(dHRUM)
library(data.table)
numdata =365
probwet =0.9
meanifwet = 8
prec= rbinom(numdata,1,probwet)*rexp(numdata,1/meanifwet)

prec[1:2] = 0.0
temp=rnorm(numdata,-4,3)
# plot(prec, type="l")
# plot(temp, type="l")
#nHrus <- 15000

# Eliades model

nHrus <- 1
ntreads <-1
#Areas <- runif(nHrus,min = 1,max  = 10) #[m2]
Areas <- runif(nHrus,min = 38780000,max  = 38780050)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus,ntreads)


setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep("DDF",times= length(Areas)),hruIds=IdsHrus)
# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus)
# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Eliades",times= length(Areas)),hruIds=IdsHrus)


setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
#setPondToAlldHrus(dHRUM_ptr = dhrus,PondTypes=rep("Pond",times= length(Areas)),hruIds=IdsHrus)
# pondDF1 = data.frame( PondArea = 40500, PonsMax= 45000, MRF= 0.039, Coflw=0.3)
# pondDF2 = data.frame( Pond_ET = "ETpond1", Pond_inSOIS= "noPondSOISPerc", Pond_inGW = "noPondGWPerc",
#                       Pond_outSOIS= "noPondSOISPerc", Pond_outGW= "PondGWPerc1",Pond_outReg="PondRouT3" )
# setPondToOnedHru(dHRUM_ptr = dhrus,0,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
# setPondToOnedHru(dHRUM_ptr = dhrus,5,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
LAdt = as.data.table(readRDS("LAI_POH_Category.rds"))
ctgrs = unique(LAdt$Category)
ctgrs
lai = LAdt[Category %in% ctgrs[17],mean_LAI]

plot(lai)

# setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
setPTLInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, Lai = lai, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 35, B_EVAP = 2.5,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.5,
                    SDIV = 0.2, CAN_ST = 2, STEM_ST = 1, CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = -2.0,
                    RETCAP = 10, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
                    KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3,CMIN =25,L=0.1, B_EXP = 0.3, KFR = 0.03,
                    INTstMax = 2,INTstScale = 1, CSfrac = 1)


setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))


# getCurdHRUpars(dHRUM_ptr = dhrus,0)
# getAllHRUpars(dHRUM_ptr = dhrus)

# getCurSHRUconfig(dHRUM_ptr = dhrus,0)
# getAllHRUconfigs(dHRUM_ptr = dhrus)



calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
# calcHBInAlldHrus(dHRUM_ptr = dhrus)
# gatherHBdata(dHRUM_ptr = dhrus)
# outDt <- getOutputDist(dHRUM_ptr = dhrus)
# dd=as.data.table(outDt$outDta)

outDta <- dHRUMrun(dHRUM_ptr = dhrus)


outDF <- data.frame(outDta$outDta)
names(outDF) <-c(outDta$VarsNams)
outDF = as.data.table(outDF)

plot(outDF$INTS, type ="l",col="green")
plot(outDF$TEMP, type ="l")
plot(outDF$PREC, type ="l")
plot(outDF$SNOW, type ="l")
plot(outDF$MELT, type ="l")

plot(outDF$PET, col ="red",t="l")
lines(outDF$EVAC, type ="l")
lines(outDF$PET, col ="red")
lines(outDF$AET, col ="blue")
lines(outDF$PET, col ="red")
lines(outDF$INTS, type ="l",col="green")

# van Dijk model
numdata =365
probwet =0.9
meanifwet = 8
prec= rbinom(numdata,1,probwet)*rexp(numdata,1/meanifwet)

prec[1:2] = 0.0
temp=rnorm(numdata,-5,3)
# plot(prec, type="l")
# plot(temp, type="l")
#nHrus <- 15000
nHrus <- 1
ntreads <-1
#Areas <- runif(nHrus,min = 1,max  = 10) #[m2]
Areas <- runif(nHrus,min = 38780000,max  = 38780050)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus,ntreads)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus)
# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Eliades",times= length(Areas)),hruIds=IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
#setPondToAlldHrus(dHRUM_ptr = dhrus,PondTypes=rep("Pond",times= length(Areas)),hruIds=IdsHrus)
# pondDF1 = data.frame( PondArea = 40500, PonsMax= 45000, MRF= 0.039, Coflw=0.3)
# pondDF2 = data.frame( Pond_ET = "ETpond1", Pond_inSOIS= "noPondSOISPerc", Pond_inGW = "noPondGWPerc",
#                       Pond_outSOIS= "noPondSOISPerc", Pond_outGW= "PondGWPerc1",Pond_outReg="PondRouT3" )
# setPondToOnedHru(dHRUM_ptr = dhrus,0,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
# setPondToOnedHru(dHRUM_ptr = dhrus,5,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))


setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 35, B_EVAP = 2.5,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.5,
                    SDIV = 0.2, CAN_ST = 2, STEM_ST = 1, CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = -2.0,
                    RETCAP = 10, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
                    KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3,CMIN =25,L=0.1, B_EXP = 0.3, KFR = 0.03,
                    INTstMax = 4,INTstScale = 1, CSfrac = 1)


setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))

calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
# calcHBInAlldHrus(dHRUM_ptr = dhrus)
# gatherHBdata(dHRUM_ptr = dhrus)
# outDt <- getOutputDist(dHRUM_ptr = dhrus)
# dd=as.data.table(outDt$outDta)

outDta <- dHRUMrun(dHRUM_ptr = dhrus)


outDF <- data.frame(outDta$outDta)
names(outDF) <-c(outDta$VarsNams)
outDF = as.data.table(outDF)

plot(outDF$INTS, type ="l")
plot(outDF$TEMP, type ="l")
plot(outDF$PREC, type ="l")

plot(outDF$PET, col ="red",t="l")
lines(outDF$EVAC, type ="l")
lines(outDF$PET, col ="red")
