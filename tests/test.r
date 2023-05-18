library(profvis)
profvis({
library(dHRUM)
numdata =36500
probwet =0.6
meanifwet = 0.1
prec= rbinom(numdata,1,probwet)*rexp(numdata,1/meanifwet)
temp=rnorm(numdata,4,3)
nHrus <- 15000
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM2",times= length(Areas)),hruIds=IdsHrus)
setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.5, KF = 0.5, ADIV = 0.5, CDIV = 0.03,
SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 0 )
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"THORNTHWAITE")
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
outDt <- getOutputDist(dHRUM_ptr = dhrus)
# outDF <- cbind(outDt$outDta, outDt$Ids)
# outDF <- data.frame(outDF)
# names(outDF) <-c(outDt$VarsNams,"HruIDs")
# outDF
})
