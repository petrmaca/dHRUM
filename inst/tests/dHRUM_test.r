sessionInfo()
library(dHRUM)
nHrus <-  200
Areas <- runif(nHrus,min = 1,max  = 100)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)

filname2 = "../dHRUM/inst/tests/indata/BP_1960_01_01.txt"
setPTInputsToAlldHrusFromFile(dHRUM_ptr = dhrus, filname2)
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Hamon")

ParDF = data.frame( B_SOIL = 1.6, C_MAX = 500, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
RETCAP = 1 )
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
for( i in 1:1){
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
}
dta <- getOutput(dHRUM_ptr = dhrus)
dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams
dF
plot(dF$TROF,type='l')
plot(dF$SURS,type='l')
plot(dF$TOTR,type='l')
plot(dF$DIRR,type='l')
plot(dF$BASF,type='l')
plot(dF$EVBS,type='l')
for(i in 1:ncol(dta$outDta)){
  print(c(min(dta$outDta[,i]),dta$VarsNams[i]))
}
plot(dF$SOIS,type='l')
plot(dF$SURS,type='l')
plot(dF$TROF,type='l')
plot(dF$GROS,type='l')

which.min(dF$SURS)
plot(dF$SURS[460:480],type='l')

nHrus <- 20
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus[1])
dhrus <- initdHruModel(nHrus,Areas[1:2],IdsHrus)
dhrus <- initdHruModel(nHrus-1,Areas[1:2],IdsHrus)

dhrus <- initdHruModel(nHrus,Areas,IdsHrus)

filname2 = "../dHRUM/inst/tests/indata/BP_1960_01_01.txt"

setInputsToAlldHrus(dHRUM_ptr = dhrus, filname2)

ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )

setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))

wrongParDF = data.frame(
  B_SOIL = 1.6,
  C_MAX = 100,
  B_EVAP = 2,
  KS = 0.1,
  KF = 0.2,
  ADIV = 0.3,
  CDIV = 0.03,
  SDIV = 0.03,
  CAN_ST = 2,
  STEM_ST = 2,
  CSDIV = 0.3,
  TETR = 5,
  DDFA = 0.5,
  TMEL = 0,
  RETCAP = 10,
  RETCAP1 = 10,
  RETCAP1 = 10
)
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(wrongParDF[1,]),names(wrongParDF))

wrongParDF1 = data.frame(
  B_SOIL = 1.6,
  C_MAX = 100,
  B_EVAP = 2,
  KS = 0.1,
  KF = 0.2,
  ADIV = 0.3,
  CDIV = 0.03,
  SDIV = 0.03,
  CAN_ST = 2,
  STEM_ST = 2,
  CSDIV = 0.3,
  TETR = 5,
  DDFA = 0.5,
  TMEL = 0,
  RETCAP1 = 10
)

setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(wrongParDF1[1,]),names(wrongParDF1))


calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Hamon")

calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Oudin")
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"erjtjkewl")
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"OudIn")

calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
outfilet=("../dHRUM/inst/tests/outdata/Rpck_dhruout.txt")
printToFile(dHRUM_ptr = dhrus,namOutFilet = outfilet)

dta<-getOutput(dHRUM_ptr = dhrus)

nHrus <- 2
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
prec=c(1,2,3)
temp=c(1,2,3)
setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, inDate = as.Date("1990/01/30"))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Oudin")
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 50, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.03, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 1 )
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
dta <- getOutput(dHRUM_ptr = dhrus)
dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams
dF
plot(dF$TROF,type='l')
plot(dF$SURS,type='l')
plot(dF$TOTR,type='l')
plot(dF$DIRR,type='l')
plot(dF$BASF,type='l')
for(i in 1:ncol(dta$outDta)){
  print(c(min(dta$outDta[,i]),dta$VarsNams[i]))
}
plot(dF$SOIS,type='l')
plot(dF$SURS,type='l')
plot(dF$TROF,type='l')
plot(dF$GROS,type='l')

# memory check
library(pryr)
mem_used()

nHrus <- 200
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
prec=c(1,2,3)
temp=c(1,2,3)
setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Oudin")
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 50, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.03, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 1 )
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
dta <- getOutput(dHRUM_ptr = dhrus)
dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams
dF
plot(dF$TROF,type='l')
plot(dF$SURS,type='l')
plot(dF$TOTR,type='l')
plot(dF$DIRR,type='l')
plot(dF$BASF,type='l')
for(i in 1:ncol(dta$outDta)){
  print(c(min(dta$outDta[,i]),dta$VarsNams[i]))
}
plot(dF$SOIS,type='l')
plot(dF$SURS,type='l')
plot(dF$TROF,type='l')
plot(dF$GROS,type='l')


library(dHRUM)
nHrus <- 10
Areas <- runif(nHrus,min = 1,max  = 100)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 500, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 1 )
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
parsmat<-matrix(1,nrow=40,ncol=40)
setParsToDistdHRUM(dHRUM_ptr = dhrus, ParsDF = data.frame(parsmat),FALSE)
parsmat<-matrix(1,nrow=40,ncol=15)
setParsToDistdHRUM(dHRUM_ptr = dhrus, ParsDF = data.frame(parsmat),TRUE)

ParDF = data.frame( C_MAX = 500, B_SOIL = 1.6, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 1 )
parsvec=as.numeric(ParDF[1,])
for(i in 1:(nHrus-1)){
  ParDF <- rbind(ParDF,parsvec)
}

ParDF[,1] <- seq(1:nHrus)
setParsToDistdHRUM(dhrus, ParDF, TRUE)
ParDF[,2] <- seq(1:nHrus)
setParsToDistdHRUM(dhrus, ParDF, FALSE)
ParDF[,3] <- seq(1:nHrus)
setParsToDistdHRUM(dhrus, ParDF, TRUE)

nHrus <- 10
Areas <- runif(nHrus,min = 1,max  = 100)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)

datDF=data.frame(HruId=c(1,2,3))
setPTInputsToDistdHRUM(dhrus,DataDF = datDF)

nHrus <- 3
Areas <- runif(nHrus,min = 1,max  = 100)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)

datDF=data.frame(HruId=c(1,2,3,3))
setPTInputsToDistdHRUM(dhrus,DataDF = datDF)

datDF=data.frame(HruId=c(1,1,2,2,3,3))
setPTInputsToDistdHRUM(dhrus,DataDF = datDF)


require(dHRUM)
nHrus <- 3
Areas <- runif(nHrus,min = 1,max  = 100)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)

datDF=data.frame(DTM=c(as.Date("2020-01-01"),
                       as.Date("2020-01-01"),
                       as.Date("2020-01-01")),
                       T=c(1,2,3),
                       P=c(20,10,1),
                       HruId=c(1,2,3))

setPTInputsToDistdHRUM(dhrus,datDF)

