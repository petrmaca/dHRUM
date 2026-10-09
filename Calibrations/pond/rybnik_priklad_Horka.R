library(dHRUM)
library(data.table)

#PRO TVORBU GRAFU
library(dygraphs)
library(xts) 
library(htmlwidgets)


source("00_functions_00.R")

################## vstupni data #######################################################################
# informace o povodich (LAT, LON,PLOCHA)
INF <- as.data.table(readRDS("./input_data/INF.rds"))#
# vstupni casove rady (teploty, srazky, LAI)
INP <- as.data.table(readRDS('./input_data/INP.rds'))#
#######################################################################################################


LIST_POV=INF$chp_14_s
SEL_POV=LIST_POV[1]

#kalibracni info (parametry + struktura) pro zvolene povodi
cal_dta=vyber_model(sel_pov=SEL_POV,hlavni_slozka = "./input_data/data_od_martina")
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


############# Model bez rybniku #############
nHrus = 1
Areas = SEL_POV_inf$a_km2*1000000 # areas je v m2
IdsHrus = paste0("MODEL_hru")
dhrus = initdHruModel(nHrus, Areas, IdsHrus,1)
setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep(SEL_SNOW,times = length(Areas)),hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep(SEL_INTERC,times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep(SEL_SmaxLAI,times= length(Areas)))
setGWtypeToAlldHrus(dhrus,gwTypes = rep(SEL_PODZ, times =nHrus), IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep(SEL_SURSTOR,times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep(SEL_FASTRESP,times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
setNumFastResAlldHrus(dHRUM_ptr = dhrus, numFastRes = POC_N, hruIds = IdsHrus)
setPTLInputsToAlldHrus(dHRUM_ptr = dhrus, Prec = SEL_POV_dta[, PREP], Temp = SEL_POV_dta[, TAVG], Lai = SEL_POV_dta[, LAI], inDate = DTM_POC)
calcPetToAllHrus(dHRUM_ptr = dhrus, Latitude = LAT, PetTypeStr = rep(SEL_PET,times = length(Areas)))

setParamsToAlldHrus(dHRUM_ptr = dhrus, as.numeric(par_values),par_names )
outDta_sim = dHRUMrun(dHRUM_ptr = dhrus)
outDF_sim = data.frame(outDta_sim$outDta)
names(outDF_sim) = outDta_sim$VarsNams



############# Model s rybnikem #############
nHrus = 1
Areas = SEL_POV_inf$a_km2*1000000
IdsHrus = paste0("MODEL_hru_pond")
dhrus = initdHruModel(nHrus, Areas, IdsHrus,1)
setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep(SEL_SNOW,times = length(Areas)),hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep(SEL_INTERC,times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep(SEL_SmaxLAI,times= length(Areas)))
setGWtypeToAlldHrus(dhrus,gwTypes = rep(SEL_PODZ, times =nHrus), IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep(SEL_SURSTOR,times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep(SEL_FASTRESP,times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
setNumFastResAlldHrus(dHRUM_ptr = dhrus, numFastRes = POC_N, hruIds = IdsHrus)
setPTLInputsToAlldHrus(dHRUM_ptr = dhrus, Prec = SEL_POV_dta[, PREP], Temp = SEL_POV_dta[, TAVG], Lai = SEL_POV_dta[, LAI], inDate = DTM_POC)
calcPetToAllHrus(dHRUM_ptr = dhrus, Latitude = LAT, PetTypeStr = rep(SEL_PET,times = length(Areas)))
setParamsToAlldHrus(dHRUM_ptr = dhrus, as.numeric(par_values),par_names )
################################ zadani rybniku #############################################################################
# pondDF1 and pondDF2 are pond variables and specifications
# PondArea = plocha nadrze
# PonsMax = maximalni objem nadrze
# MRF = minimalni zustatkovy prutok
# Coflw = konstantní odtok
# Pond_ET = "ETpond1" = volba rovnice pro vypar z hladiny - A. Beran
# Pond_outReg="PondRouT3" =  Coflw - konstantni odtok při této volbě přebírá hodnotu Coflw
# Pond_inSOIS,Pond_inGW, Pond_outSOIS,Pond_outGW - prusaky nadrzi/nyní vypnuty

PondArea_m2 = 45000
PondVolume_m3 =45000
MRF= 0 
Coflw=0.03 # pokud dam coflw = 0 tak nejde vydet kolisani objemu nadrze. nadrz se naplni a je porad plna


pondDF1 = data.frame( PondArea = PondArea_m2, PonsMax= PondVolume_m3, MRF= MRF, Coflw=Coflw)
pondDF2 = data.frame( Pond_ET = "ETpond1", Pond_inSOIS= "noPondSOISPerc", Pond_inGW = "noPondGWPerc",
                      Pond_outSOIS= "noPondSOISPerc", Pond_outGW= "noPondGWPerc",Pond_outReg="PondRouT3" )
# pond implementation
setPondToOnedHru(dHRUM_ptr = dhrus,0,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
################################ test pond - end #############################################################################
outDta_sim_pond = dHRUMrun(dHRUM_ptr = dhrus)
outDF_sim_pond = data.frame(outDta_sim_pond$outDta)
names(outDF_sim_pond) = outDta_sim_pond$VarsNams


################################ srovnavaci grafy #############################################################################

graf_porovnani_TOTR(outDF_sim_pond,PondArea_m2,CatchmentArea_m2=Areas,uloz_graf=TRUE)
graf_porovnani_pritoku_zasoby(outDF_sim_pond,PondArea_m2,CatchmentArea_m2=Areas,uloz_graf=TRUE)
graf_vypar_z_nadrze(outDF_sim_pond,PondArea_m2,uloz_graf=TRUE)
graf_porovnani_vyparu(outDF_sim,outDF_sim_pond,PondArea_m2,CatchmentArea_m2=Area,uloz_graf=TRUE)









