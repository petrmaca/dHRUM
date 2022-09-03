library(dHRUM)
nHrus <- 1
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
prec=c(100,200,300,100,200,300)
temp=c(1,2,3,1,2,3)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.5, KF = 0.5, ADIV = 0.5, CDIV = 0.03,
SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 0 )
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"THORNTHWAITE")
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
outDta <- getOutputDist(dHRUM_ptr = dhrus)
outDF <- cbind(outDta$outDta, outDta$Ids)
outDF <- data.frame(outDF)
names(outDF) <-c(outDta$VarsNams,"HruIDs")
outDF

