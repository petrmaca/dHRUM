# library(profvis)
# profvis(
# {
library(dHRUM)
library(data.table)
numdata =36500
probwet =0.9
meanifwet = 8
prec= rbinom(numdata,1,probwet)*rexp(numdata,1/meanifwet)
temp=rnorm(numdata,20,3)
#nHrus <- 15000
nHrus <- 1
#Areas <- runif(nHrus,min = 1,max  = 10) #[m2]
Areas <- runif(nHrus,min = 38780000,max  = 38780050)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus,64)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)

#setPondToAlldHrus(dHRUM_ptr = dhrus,PondTypes=rep("Pond",times= length(Areas)),hruIds=IdsHrus)
pondDF1 = data.frame( PondArea = 40500, PonsMax= 45000, MRF= 0.039)
pondDF2 = data.frame( Pond_ET = "ETpond1", Pond_inSOIS= "noPondSOISPerc", Pond_inGW = "noPondGWPerc",
                      Pond_outSOIS= "noPondSOISPerc", Pond_outGW= "PondGWPerc3",Pond_outReg="PondRouT1" )

setPondToOnedHru(dHRUM_ptr = dhrus,0,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))

setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 35, B_EVAP = 2.5,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.2,
                    SDIV = 0.1, CAN_ST = 2, STEM_ST = 1, CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 10, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
                    KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3,CMIN =25,L=0.1, B_EXP = 0.3, KFR = 0.03)
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))

current_parameters<-getCurdHRUpars(dHRUM_ptr = dhrus,0)

calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
# calcHBInAlldHrus(dHRUM_ptr = dhrus)
# gatherHBdata(dHRUM_ptr = dhrus)
# outDt <- getOutputDist(dHRUM_ptr = dhrus)
# dd=as.data.table(outDt$outDta)

outDta <- dHRUMrun(dHRUM_ptr = dhrus)


outDF <- data.frame(outDta$outDta)
names(outDF) <-c(outDta$VarsNams)
outDF = as.data.table(outDF)
#outDF$SOIS
#outDF$AET
#outDF$PET
#outDF$EVBS
#outDF$EVAS
#outDF$EVAC
#outDF$ETWS
#outDF$PONS
#tail(outDF$TOTR)


max(outDF$AET - outDF$PET)
# )
smax =(ParDF$B_SOIL*ParDF$CMIN + ParDF$C_MAX) / (ParDF$B_SOIL+1)
plot(outDF$SOIS, type ="l")
abline(h=smax,col="red")
plot(outDF$SURS, type ="l")
plot(outDF$DIRR, type ="l")
plot(outDF$PONS, type ="l")
plot(outDF$TOTR, type ="l")
