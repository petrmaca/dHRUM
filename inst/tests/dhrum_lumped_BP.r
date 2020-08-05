library(dHRUM)
nHrus <- 1
Areas <- 4300000
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhruBP <- initdHruModel(nHrus,Areas,IdsHrus)

filname2 = "../dHRUM/inst/tests/indata/BP_1960_01_01.txt"
setInputsToAlldHrus(dhruBP, filname2)
# Hamon PET
calcPetToAllHrus(dhruBP,50.1,"Hamon")
# Oudin PET
# calPetToAllHrus(dHRU_ptr = dhruBP,50.1,"Oudin")

attach(what = "data/Amalie_lumped_dHRUM.rda")
Par_dHRUm_Amalie_lumped$pars[1,]
ParBP = Par_dHRUm_Amalie_lumped$pars[1,1:15]
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

