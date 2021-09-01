library(RcppDE)
library(data.table)
library(dHRUM)

NhrusBP <- nrow(br)

dny=c(30,60,90,120,150,180,210,240,270,300,330,355,364)
p_OBS=dny/365.25

# BP 96.68961 (1.736-2.430446)
RaBP = 96# odhad Martin Hanel
QmBP = c(26, 18, 14, 12, 10, 8.0, 7.0, 6.0, 4.5, 3.5, 2.5, 1.0, 0.5)
A=sum(SoilBP$Area)# plocha BP
RmBP = QmBP * (3600*24) / A #CHMU ZHU 4.42 v datech SoilBP 4.41

nHrus <- NhrusBP
Areas <- br$Plocha
IdsHrus <- br$ID
dhrusBP <- initdHruModel(nHrus,Areas,IdsHrus)
filname2 = "../dHRUM/Calibrations/Amalie/indata/BP_1960_01_01.txt"
setPTInputsToAlldHrusFromFile(dHRUM_ptr = dhrusBP, filname2)
calcPetToAllHrus(dHRUM_ptr = dhrusBP,50.1,"HAMON")
ParDFup1 = data.frame( B_SOIL = 2., C_MAX = 800, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.1,
                       SDIV = 0.05, CAN_ST = 2., STEM_ST = 2., CSDIV = 0.8, TETR = 4, DDFA = 0.8, TMEL = 0.0,
                       RETCAP = 2 )

parsvec=as.numeric(ParDFup1[1,])
for(i in 1:(nHrus-1)){
  ParDFup1 <- rbind(ParDFup1,parsvec)
}
ParDFlow1 = data.frame( B_SOIL = 0.6, C_MAX = 10, B_EVAP = 0.6,  KS = 0.002, KF = 0.2, ADIV = 0.01, CDIV = 0.05,
                        SDIV = 0.01, CAN_ST = 0.1, STEM_ST = 0.1, CSDIV = 0.01, TETR = 0, DDFA = 0.08, TMEL = -10.0,
                        RETCAP = 1 )

parsvec=as.numeric(ParDFlow1[1,])
for(i in 1:(nHrus-1)){
  ParDFlow1 <- rbind(ParDFlow1,parsvec)
}
ParNams <- names(ParDFlow1)

nParIhru <- ncol(ParDFlow1)

mae = function(myPar){
  newmat <- as.data.frame(matrix(myPar,nrow = nHrus, ncol = nParIhru))
  names(newmat) <- ParNams
  setParsToDistdHRUM(dhrusBP, newmat, FALSE)
  # # for( i in 1:1000){
  calcHBInAlldHrus(dHRUM_ptr = dhrusBP)
  gatherHBdata(dHRUM_ptr = dhrusBP)
  # # }
  dta <- getOutput(dHRUM_ptr = dhrusBP)
  dF <-data.frame(dta$outDta)
  names(dF) <- dta$VarsNams
  simRM=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
  mymae =as.double(sum(abs(simRM - RmBP)))
  # # mymae=NA
  # if(is.na(mymae)) mymae = 9999
  # (return as.double(mymae))
}

itermaxW=2
decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 100,
                        itermax = itermaxW, CR = 0.25, F = 0.7, trace = TRUE,
                        initialpop = NULL, storepopfrom = itermaxW + 1,
                        storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                        steptol = itermaxW)

u=DEoptim( lower=as.numeric(as.matrix(ParDFlow1)),
           upper=as.numeric(as.matrix(ParDFup1)), fn=mae, control = decntr)

# u$optim$bestmem

ParBestVec <- as.numeric(u$optim$bestmem)
ParBestDF <- as.data.frame(matrix(ParBestVec,nrow = nHrus, ncol = nParIhru))
names(ParBestDF) <- ParNams
setParsToDistdHRUM(dhrusBP, ParBestDF, TRUE)
dta <- dHRUMrun(dHRUM_ptr = dhrusBP)

dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams
plot(dF$TOTR, type='l')
plot(dF$BASF, type='l')
plot(dF$DIRR, type='l')
plot(dF$SOIS, type='l')
plot(dF$GROS, type='l')
plot(dF$SNOW, type='l')
simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
plot(RmBP, simBest)

plot(p_OBS,RmBP, pch=19, ylim=range(c(RmBP,simBest)),
     main = paste0("Brejlský potok  - ",nHrus," HRU jednotek"),xlab="P(Qm)", ylab="Qm [mm/den]")
points(p_OBS,simBest,col="red",pch=19)
grid()
mm=mae(as.numeric(u$optim$bestmem))
mm


NhrusKL <- nrow(kl)


dny=c(30,60,90,120,150,180,210,240,270,300,330,355,364)
p_OBS=dny/365.25

RaKL = 96# odhad Martin Hanel
# ZHU Qa
Qa=10
# mdenni Q
QmKL = c(22, 15, 12, 10, 8.5, 6.5, 6.0, 5.0, 3.5, 3.0, 2.0, 1.0, 0.5)
A=sum(SoilKL$Area)# plocha KL 3.28 chmu plocha
RmKL = QmKL * (3600*24) / A #CHMU ZHU
Ra =Qa*3600*24*365/A

RmKL = QmKL * (3600*24) / A #CHMU ZHU 3.28 v datech SoilKL 3.24


filname2 = "../dHRUM/Calibrations/Amalie/indata/KL_1960_01_01_noDate.txt"
Areas <- kl$Plocha
IdsHrus <- kl$ID
dhrusKL <- initdHruModel(NhrusKL,Areas,IdsHrus)
PTdta = read.table(filname2)
prec=PTdta$V1
temp=PTdta$V2
setPTInputsToAlldHrus(dHRUM_ptr = dhrusKL, Prec = prec, Temp = temp, as.Date("1990/01/30"))
calcPetToAllHrus(dHRUM_ptr = dhrusKL,50.1,"HAMON")

ParDFup1 = data.frame( B_SOIL = 2., C_MAX = 800, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.1,
                       SDIV = 0.05, CAN_ST = 2., STEM_ST = 2., CSDIV = 0.8, TETR = 4, DDFA = 0.8, TMEL = 0.0,
                       RETCAP = 2 )

parsvec=as.numeric(ParDFup1[1,])
for(i in 1:(nHrus-1)){
  ParDFup1 <- rbind(ParDFup1,parsvec)
}
# ParDFup1$C_MAX <- SoilBP$RVK
# ParDFup1$C_MAX <- SoilBP$VVK

ParDFlow1 = data.frame( B_SOIL = 0.6, C_MAX = 10, B_EVAP = 0.6,  KS = 0.002, KF = 0.2, ADIV = 0.01, CDIV = 0.05,
                        SDIV = 0.01, CAN_ST = 0.1, STEM_ST = 0.1, CSDIV = 0.01, TETR = 0, DDFA = 0.08, TMEL = -10.0,
                        RETCAP = 1 )

parsvec=as.numeric(ParDFlow1[1,])
for(i in 1:(nHrus-1)){
  ParDFlow1 <- rbind(ParDFlow1,parsvec)
}


ParNams <- names(ParDFlow1)

nParIhru <- ncol(ParDFlow1)
mae = function(myPar){
  newmat <- as.data.frame(matrix(myPar,nrow = nHrus, ncol = nParIhru))
  names(newmat) <- ParNams
  setParsToDistdHRUM(dhrusKL, newmat, FALSE)
  # # for( i in 1:1000){
  calcHBInAlldHrus(dHRUM_ptr = dhrusKL)
  gatherHBdata(dHRUM_ptr = dhrusKL)
  # # }
  dta <- getOutput(dHRUM_ptr = dhrusKL)
  dF <-data.frame(dta$outDta)
  names(dF) <- dta$VarsNams
  simRM=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
  mymae =as.double(sum(abs(simRM - RmKL)))
  # # mymae=NA
  # if(is.na(mymae)) mymae = 9999
  # (return as.double(mymae))
}

library(RcppDE)
itermaxW=2
decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 200,
                        itermax = itermaxW, CR = 0.25, F = 0.7, trace = TRUE,
                        initialpop = NULL, storepopfrom = itermaxW + 1,
                        storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                        steptol = itermaxW)

u=DEoptim( lower=as.numeric(as.matrix(ParDFlow1)),
           upper=as.numeric(as.matrix(ParDFup1)), fn=mae, control = decntr)

# u$optim$bestmem

ParBestVec <- as.numeric(u$optim$bestmem)
ParBestDF <- as.data.frame(matrix(ParBestVec,nrow = nHrus, ncol = nParIhru))
names(ParBestDF) <- ParNams

setParsToDistdHRUM(dhrusKL, ParBestDF, TRUE)
dta <- dHRUMrun(dHRUM_ptr = dhrusKL)

dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams
plot(dF$TOTR, type='l')
plot(dF$BASF, type='l')
plot(dF$DIRR, type='l')
plot(dF$SOIS, type='l')
plot(dF$GROS, type='l')
plot(dF$SNOW, type='l')
simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
plot(RmBP, simBest)
simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
plot(RmKL, simBest)

plot(p_OBS,RmKL, pch=19, ylim=range(c(RmKL,simBest)),
     main = paste0("Karlův Luh - ",NhrusKL," HRU jednotek"),xlab="P(Qm)", ylab="Qm [mm/den]")
points(p_OBS,simBest,col="red",pch=19)
grid()

mmKL=mae(as.numeric(u$optim$bestmem))
mmKL
