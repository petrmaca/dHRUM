library(dHRUM)
library(data.table)

#PRO TVORBU GRAFU
library(dygraphs)
library(xts)
library(htmlwidgets)

source("C:/Users/Lenovo/Documents/dHRUM/Calibrations/wetland/00_functions_00.R")

################## vstupni data #######################################################################
# informace o povodich (LAT, LON,PLOCHA)
INF <- as.data.table(readRDS("C:/Users/Lenovo/Documents/dHRUM/Calibrations/pond/input_data/INF.rds"))#
# vstupni casove rady (teploty, srazky, LAI)
INP <- as.data.table(readRDS('C:/Users/Lenovo/Documents/dHRUM/Calibrations/pond/input_data/INP.rds'))#
#######################################################################################################INF


LIST_POV=INF$chp_14_s
SEL_POV=LIST_POV[1]

#kalibracni info (parametry + struktura) pro zvolene povodi
cal_dta=vyber_model(sel_pov=SEL_POV,hlavni_slozka = "C:/Users/Lenovo/Documents/dHRUM/Calibrations/pond/input_data/data_od_martina")
parametry=cal_dta[,25:36]
par_names=names(parametry)
par_values=as.numeric(parametry)

#omezení vstupnich dat na vybrane povodi
SEL_POV_dta = INP[INP$ID %in% SEL_POV]
SEL_POV_inf = INF[INF$chp_14_s %in% SEL_POV]

DTM_POC = SEL_POV_dta[1, DTM]
LAT = SEL_POV_inf$LAT

# nastaveni modelu pro dane povodi / definice struktury modelu
SEL_SNOW = cal_dta$SNIH
SEL_PET = cal_dta$PET
SEL_INTERC = cal_dta$INTER
SEL_SmaxLAI = cal_dta$MET_INTR_MAX
SEL_PODZ = cal_dta$PODZ
SEL_SURSTOR = cal_dta$POV_RET
SEL_FASTRESP = cal_dta$FRS
SEL_PUDA = cal_dta$PUDA
POC_N = cal_dta$NUM_FSR


############# Model bez wtlnd #############
nHrus = 1
Areas = SEL_POV_inf$a_km2*1000000 # areas je v m2
IdsHrus = paste0("MODEL_hru")
dhrus = initdHruModel(nHrus, Areas, IdsHrus,1)
setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep(SEL_SNOW,times = length(Areas)),hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep(SEL_INTERC,times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep(SEL_SmaxLAI,times= length(Areas)))
setGWtypeToAlldHrus(dhrus,gwTypes = rep(SEL_PODZ, times =nHrus), IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep(SEL_SURSTOR,times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep(SEL_FASTRESP,times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep(SEL_PUDA,times= length(Areas)), hruIds = IdsHrus)
setNumFastResAlldHrus(dHRUM_ptr = dhrus, numFastRes = POC_N, hruIds = IdsHrus)
setPTLInputsToAlldHrus(dHRUM_ptr = dhrus, Prec = SEL_POV_dta[, PREP], Temp = SEL_POV_dta[, TAVG], Lai = SEL_POV_dta[, LAI], inDate = DTM_POC)
calcPetToAllHrus(dHRUM_ptr = dhrus, Latitude = LAT, PetTypeStr = rep(SEL_PET,times = length(Areas)))

setParamsToAlldHrus(dHRUM_ptr = dhrus, as.numeric(par_values),par_names )
outDta_sim = dHRUMrun(dHRUM_ptr = dhrus)
outDF_sim = data.frame(outDta_sim$outDta)
names(outDF_sim) = outDta_sim$VarsNams



############# Model s wtlnd #############

nHrus = 1
Areas = SEL_POV_inf$a_km2*1000000
IdsHrus = paste0("MODEL_hru_pond")
dhrus = initdHruModel(nHrus, Areas, IdsHrus,1)

Wetland_area=Areas*0.01  # Wetland area is 1% of catchment area
Wpar_names=c("WtlnFrac","Kinct","KwPe")
Wpar_values=c(Wetland_area/Areas,0.1,0.2)
SEL_SURSTOR ="Wetland"
# WtlnFrac: relative value of wetland size [0,1]
# Kinct: Beer-Lambert extinction coefficient [0,1]
# KwPe: Storage coefficient for wetland percolation [0,1]

wpar_values=c(par_values,Wpar_values)
wpar_names=c(par_names,Wpar_names)

setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep(SEL_SNOW,times = length(Areas)),hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep(SEL_INTERC,times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep(SEL_SmaxLAI,times= length(Areas)))
setGWtypeToAlldHrus(dhrus,gwTypes = rep(SEL_PODZ, times =nHrus), IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep(SEL_SURSTOR,times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep(SEL_FASTRESP,times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep(SEL_PUDA,times= length(Areas)), hruIds = IdsHrus)
setNumFastResAlldHrus(dHRUM_ptr = dhrus, numFastRes = POC_N, hruIds = IdsHrus)
setPTLInputsToAlldHrus(dHRUM_ptr = dhrus, Prec = SEL_POV_dta[, PREP], Temp = SEL_POV_dta[, TAVG], Lai = SEL_POV_dta[, LAI], inDate = DTM_POC)
calcPetToAllHrus(dHRUM_ptr = dhrus, Latitude = LAT, PetTypeStr = rep(SEL_PET,times = length(Areas)))

setParamsToAlldHrus(dHRUM_ptr = dhrus, as.numeric(wpar_values),wpar_names )

outDta_sim_measure = dHRUMrun(dHRUM_ptr = dhrus)
outDF_sim_measure= data.frame(outDta_sim_measure$outDta)
names(outDF_sim_measure) = outDta_sim_measure$VarsNams



################################ srovnavaci grafy #############################################################################
graf_porovnani_TOTR(outDF_sim,outDF_sim_measure,CatchmentArea_m2=Areas,uloz_graf=FALSE)
graf_porovnani_celkoveho_vyparu(outDF_sim,outDF_sim_measure,uloz_graf=FALSE)
graf_porovnani_vyparu_z_retence(outDF_sim,outDF_sim_measure,uloz_graf=FALSE)
graf_transpirace_WTLND(outDF_sim_measure,CatchmentArea_m2=Areas,uloz_graf=FALSE)









