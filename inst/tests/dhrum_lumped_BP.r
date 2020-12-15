library(dHRUM)
nHrus <- 1
Areas <- 4300000
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhruBP <- initdHruModel(nHrus,Areas,IdsHrus)

filname2 = "../dHRUM/inst/tests/indata/BP_1960_01_01.txt"


setPTInputsToAlldHrusFromFile(dhruBP, filname2)
# Hamon PET
calcPetToAllHrus(dhruBP,50.1,"Oudin")
# Oudin PET
# calPetToAllHrus(dHRU_ptr = dhruBP,50.1,"Oudin")

attach(what = "data/Amalie_lumped_dHRUM.rda")
Par_dHRUm_Amalie_lumped$pars[1,]
ParBP = Par_dHRUm_Amalie_lumped$pars[1,1:15]
ParBP=BP_df[1,1:15]
setParamsToAlldHrus(dhruBP,as.numeric(ParBP[1,]),names(ParBP))
calcHBInAlldHrus(dhruBP)
gatherHBdata(dhruBP)


dta <- getOutput(dhruBP)
dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams


# celkovy odtok
plot(dF$TOTR, type='l')
# zakladni odtok
plot(dF$BASF, type='l')
# primy odtok
plot(dF$DIRR, type='l')
# zasoba vody v pude
plot(dF$SOIS, type='l')
# zasoba vody v podzemni vode
plot(dF$GROS, type='l')


VecDATE = c()
for(i in 1:nrow(dF)){
  VecDATE[i] = paste0(as.character(dF[i,1]),"-",as.character(dF[i,2]),"-",as.character(dF[i,3]))
}
VecDATE=as.Date(VecDATE)

Prec=dF$PREC
Temp=dF$TEMP

nHrus <- 1
Areas <- 4300000
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhruBP1 <- initdHruModel(nHrus,Areas,IdsHrus)
setPTDateInputsToAlldHrus(dhruBP1,Prec = Prec,Temp = Temp,DateVec = VecDATE)
calcPetToAllHrus(dhruBP1,Latitude = 50.1,PetTypeStr = "Hamon")
setParamsToAlldHrus(dhruBP1,as.numeric(ParBP[1,]),names(ParBP))
calcHBInAlldHrus(dhruBP1)
gatherHBdata(dhruBP1)
dta1 <- getOutput(dhruBP1)
dF1 <-data.frame(dta1$outDta)
names(dF1) <- dta1$VarsNams

plot(dF$PREF,dF1$PREF)
sum(dF$PREF-dF1$PREF)

library(data.table)
DFDT=as.data.table(dF)

DFDT[,DTM:=as.Date(paste(YEAR,MONTH,DAY,sep="-")),]
dfShort=DFDT[DTM>"1980-01-01" & DTM<"2011-01-01",]
dtM=dfShort[,.(PRM=sum(PREC),AETM=sum(AET)),by=.(MONTH,YEAR)]

dtM[,.(median(PRM)),by=.(MONTH)]

dta=readRDS("./inst/tests/indata/StepProSimulace.rds")
dtaBP=dta[DTM>"1979-12-31" &  DTM<"2011-01-01" & ID=='BP',]

nHrus <- 1
Areas <- 4300000
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhruBP1 <- initdHruModel(nHrus,Areas,IdsHrus)
setPTDateInputsToAlldHrus(dhruBP1,Prec = dtaBP$P,Temp = dtaBP$T,DateVec = dtaBP$DTM)
calcPetToAllHrus(dhruBP1,Latitude = 50.1,PetTypeStr = "Oudin")

ParBest[1,] = as.numeric(BP_df[1,1:15])
ParBest

BP_df


setParamsToAlldHrus(dhruBP1,as.numeric(ParBest[1,]),names(ParBest))
calcHBInAlldHrus(dhruBP1)
gatherHBdata(dhruBP1)
dta1 <- getOutput(dhruBP1)
dF <-data.frame(dta1$outDta)
names(dF) <- dta1$VarsNams
DFDT=as.data.table(dF)

DFDT[,DTM:=as.Date(paste(YEAR,MONTH,DAY,sep="-")),]
DFDT[,DTM:=as.Date(paste(YEAR,MONTH,DAY,sep="-")),]
dfShort=DFDT[DTM>"1980-01-01" & DTM<"2011-01-01",]
# na.omit(dfShort)
dfShort=DFDT[DTM>"1980-01-01" & DTM<"2011-01-01",Aet:=AET+EVBS+EVAS]
dfShort=na.omit(dfShort)
dtM=dfShort[,.(PRM=sum(PREC),AETM=sum(Aet),PETM=sum(PET),SO=mean(SOIS),D=sum(DIRR),BF=sum(BASF),R=sum(TOTR)),by=.(MONTH,YEAR)]
# dtM=na.omit(dtM)
MdatBP=dtM[,.(P=median(PRM),PET=median(PETM),AET=median(AETM),SOIL=median(SO),DR=median(D),BF=median(BF),R=median(R)),by=.(MONTH)]
MdatBP
saveRDS(file = "../../bp_m.rds",MdatBP)
