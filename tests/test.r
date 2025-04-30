# library(profvis)
# profvis(
# {
library(dHRUM)
numdata =36500
probwet =0.6
meanifwet = 8
prec= rbinom(numdata,1,probwet)*rexp(numdata,1/meanifwet)
temp=rnorm(numdata,10,3)
nHrus <- 1000
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus,1)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("EXP_LOG", times =nHrus), IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM2",times= length(Areas)),hruIds=IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus)
setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 10, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
                    KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3,CMIN =20,L=0.1, B_EXP = 0.3)
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
# calcHBInAlldHrus(dHRUM_ptr = dhrus)
# gatherHBdata(dHRUM_ptr = dhrus)
outDt <- getOutputDist(dHRUM_ptr = dhrus)
library(data.table)
dd=as.data.table(outDt$outDta)

outDta <- dHRUMrun(dHRUM_ptr = dhrus)
outDF <- data.frame(outDta$outDta)
names(outDF) <-c(outDta$VarsNams)
# }
# )

