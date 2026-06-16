# library(profvis)
# profvis(
# {
library(dHRUM)
library(data.table)
numdata = 365
probwet = 0.9
meanifwet = 8
prec= rbinom(numdata,1,probwet) * rexp(numdata,1/meanifwet)

# prec[1:2] = 0.0
temp=rnorm(numdata,-1,2)
# plot(prec, type="l")
# plot(temp, type="l")
#nHrus <- 15000

# Eliades model

nHrus <- 1
ntreads <- 1
#Areas <- runif(nHrus,min = 1,max  = 10) #[m2]
Areas <- runif(nHrus,min = 38780000,max  = 38780050)
IdsHrus <- paste0("ID",seq(1:length(Areas)))


dhrus <- initdHruModel(nHrus,Areas,IdsHrus,ntreads)


setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep("DDF",times = length(Areas)),hruIds = IdsHrus)
# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times = length(Areas)),hruIds = IdsHrus,InstStLai = rep(TRUE,times = length(Areas)),smaxlaiTypes = rep("Pitman",times = length(Areas)))

# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("Pitman",times= length(Areas)))

setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Eliades",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("Pitman",times= length(Areas)))

# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("VonHoyningenHuene",times= length(Areas)))
# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("VonHoyningenHuene",times= length(Areas)))
# setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Eliades",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("VonHoyningenHuene",times= length(Areas)))


setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)

# setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
# setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfacePRTL",times= length(Areas)),hruIds=IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("Wetland",times= length(Areas)),hruIds=IdsHrus)

setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
#setPondToAlldHrus(dHRUM_ptr = dhrus,PondTypes=rep("Pond",times= length(Areas)),hruIds=IdsHrus)
# pondDF1 = data.frame( PondArea = 40500, PonsMax= 45000, MRF= 0.039, Coflw=0.3)
# pondDF2 = data.frame( Pond_ET = "ETpond1", Pond_inSOIS= "noPondSOISPerc", Pond_inGW = "noPondGWPerc",
#                       Pond_outSOIS= "noPondSOISPerc", Pond_outGW= "PondGWPerc1",Pond_outReg="PondRouT3" )
# setPondToOnedHru(dHRUM_ptr = dhrus,0,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
# setPondToOnedHru(dHRUM_ptr = dhrus,5,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
# LAdt = as.data.table(readRDS("LAI_POH_Category.rds"))
LAdt = as.data.table(get_Lai_DataCat())
ctgrs = unique(LAdt$Category)
ctgrs
lai = LAdt[Category %in% ctgrs[17],mean_LAI]

plot(lai)

# setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
setPTLInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, Lai = lai, as.Date("1990/01/30"))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 20, B_EVAP = 2.5,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.5,
                    SDIV = 0.2, CAN_ST = 2, STEM_ST = 1, CSDIV = 0.8, TETR = 0, DDFA = 6, TMEL = -1.0,
                    RETCAP = 8, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
                    KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3,CMIN =0,L=0.1, B_EXP = 0.3, KFR = 0.03,
                    INTstMax = 2,INTstScale = 1, CSfrac = 1, SRFrac =0.75,Kinct = 0.8, KwPe=0.01, Csnow = 0.5)

# (1.6*10+100)/(1.6+1) smax cmin=10 cmax =100 bsoil= 1.6
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


PrecNoSnow = outDF$PREC - outDF$SNOW

sum(PrecNoSnow) - sum(outDF$RAIN)- sum(outDF$SUBL)

sum(outDF$PREC)
sum(outDF$PREF)
sum(outDF$EVAC) + sum(outDF$PREF) + sum(outDF$SUBL) +sum(outDF$EVAS)
# sum(PrecNoSnow) + sum(outDF$SNOW) + sum(outDF$SUBL)

sum(outDF$PREC)
sum(outDF$RAIN) + sum(outDF$SNOW) + sum(outDF$SUBL)

sum(outDF$SNOW)
# sum(outDF$SNOW) - (sum(outDF$MELT) + sum(outDF$MELV) + sum(outDF$SUBL))

sum(outDF$MELT) + sum(outDF$MELV)   + sum(outDF$SUBL)
sum(outDF$REFR)

sum(outDF$TROF)

sum(PrecNoSnow) + sum(outDF$MELT) + sum(outDF$MELV) - sum(outDF$PREF)
sum(outDF$EVAC)+sum(outDF$EVAS) + sum(outDF$REFR)

sum(outDF$EVAS) + sum(outDF$EVAC)+sum(outDF$PREF)+ sum(outDF$SUBL)

sum(outDF$PREC)- (sum(outDF$PREF)+sum(outDF$EVAC)+sum(outDF$EVAS)+sum(outDF$SUBL))
# sum(outDF$SNOW-outDF$MELT)
outDF$PREF[length(outDF$PREF)]

min(outDF$PREC -outDF$PREF)
min(outDF$PREC -outDF$TROF)
min(outDF$PREF -outDF$TROF)

plot(outDF$INTS, type ="l",col="green")
plot(outDF$TEMP, type ="l")
plot(outDF$PREC, type ="l")
plot(outDF$SNOW, type ="l")
plot(outDF$MELT, type ="l")
plot(outDF$MELV, type ="l")
plot(outDF$SUBL, type ="l")

plot(outDF$INTS, type ="l",col="green",ylim=c(0,range(c(outDF$PET,outDF$INTS))[2]))
lines(outDF$EVAC, type ="l")
lines(outDF$EVAS, type ="l", col="grey")
lines(outDF$PET, col ="red")
lines(outDF$AET, col ="blue")
lines(outDF$PET, col ="red")

lines(outDF$INTS, type ="l",col="green")

plot(outDF$INTS, type ="l",col="green", ylim=range(outDF$TROF))
lines(outDF$TROF, col ="grey")
lines(outDF$PREF, col ="blue")
lines(outDF$PREC, col ="black")

lines(outDF$EVAC, type ="l")
lines(outDF$PET, col ="red")
lines(outDF$AET, col ="blue")
lines(outDF$PET, col ="red")
lines(outDF$INFL, col ="lightblue")
lines(outDF$SOIS, col ="red")
lines(outDF$PERC, col ="red")

plot(outDF$SURS, type="l")
sum(outDF$SURS)
lines(outDF$INFL,col="red")
plot(outDF$INFL,col="red", type='l')
sum(outDF$INFL)
sum(outDF$PREF)

plot(outDF$SURS,col="red",type="l", ylim=c(0, max(outDF$SURS)))
lines(outDF$TRNS, col="green")
lines(outDF$ETSW, col="blue")
lines((outDF$EVBS+outDF$TRNS), col="grey")


sum(outDF$PET)
sum(outDF$AET)
sum(outDF$ETWS)
sum(outDF$PREF)
sum(outDF$TRNS)
sum(outDF$INFL) +sum(outDF$ETWS)+sum(outDF$TRNS)

sum(outDF$PREF) - (sum(outDF$INFL) + sum(outDF$ETWS) + sum(outDF$TRNS))
outDF$SURS[length(outDF$SURS)]


sum(outDF$INFL) +sum(outDF$TRNS)+sum(outDF$ETWS)

min(outDF$PREC - outDF$PREF)
min(outDF$PREC - outDF$MELT)
min(outDF$PREC - outDF$TROF)
min(outDF$PREF - outDF$TROF)

min(outDF$PREF - outDF$INFL)

dt=data.frame(pf=outDF$PREF, tt=outDF$TEMP, sn=outDF$SNOW, m=outDF$MELT,dd=outDF$PREC -outDF$MELT)
dT=data.frame(pf=outDF$PREF, inl=outDF$INFL, e=outDF$ETWS, s=outDF$SURS,dd=outDF$PREF -outDF$INFL)

sum(outDF$TROF)

sum(outDF$PREC)
sum(outDF$PREF)
sum(outDF$EVAC)
sum(outDF$EVAS)

PrecNoSnow = outDF$PREC -outDF$SNOW

sum(PrecNoSnow) + sum(outDF$MELT)

sum(PrecNoSnow) + sum(outDF$MELT) + sum(outDF$MELV)- sum(outDF$PREF)
sum(outDF$EVAC)+sum(outDF$EVAS)

sum(outDF$PREC)- (sum(outDF$PREF)+sum(outDF$EVAC)+sum(outDF$EVAS))
sum(outDF$SNOW-outDF$MELT)
outDF$SURS[length(outDF$PREF)]

outDF$SURS[length(outDF$INTS)]

sum(outDF$SNOW)

sum(outDF$MELT)

sum(outDF$PREF)
sum(outDF$ETWS)
sum(outDF$INFL)

sum(outDF$SURS)
max(outDF$SOIS)

sum(outDF$ETSW)

sum(outDF$INFL) + sum(outDF$ETWS)

sum(outDF$PREF)

plot(outDF$SOIS, col ="red", type="l")


