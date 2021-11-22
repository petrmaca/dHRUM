library(dHRUM)
nHrus <- 1
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
prec=c(1,2,3)
temp=c(1,2,3)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, c("COLLIE_V2"), c("ID1"))
gwStorTypesVec <- c("LINBY_RES")
setGWtypeToAlldHrus(dHRUM_ptr = dhrus, gwStorTypesVec, c("ID1"))
setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
                    SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10,
                    D_BYPASS = 0.8, THR = 10, KS2 = 0.1, FOREST_FRACT = 0.3,
                    FC = 100)
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
#calcHBInAlldHrus(dHRUM_ptr = dhrus)
#gatherHBdata(dHRUM_ptr = dhrus)
outDF <- dHRUMrun(dHRUM_ptr = dhrus)
outDf <- data.frame(outDF$outDta)
names(outDf) <- outDF$VarsNams
plot(outDf$SOIS, type = 'l')
plot(outDf$EVBS, type = 'l')
plot(outDf$PERC, type = 'l')

