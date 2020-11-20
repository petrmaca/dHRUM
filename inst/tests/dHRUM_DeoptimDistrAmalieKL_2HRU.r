sessionInfo()
library(dHRUM)
nHrus <- 2
# Areas <- runif(nHrus,min = 1,max  = 100)
Areas <- c(2.23,1.05)#les 2.23 Zem_puda 1.05
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)

filname2 = "../tests/indata/KL_1960_01_01_noDate.txt"
TPdta = read.table(filname2)

prec=TPdta$V1
temp=TPdta$V2


setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1960/01/01"))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Oudin")

ParDFup1 = data.frame( B_SOIL = 2, C_MAX = 200, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.3,
                       SDIV = 0.3, CAN_ST = 4., STEM_ST = 4., CSDIV = 0.8, TETR = 0.5, DDFA = 10, TMEL = 0.0,
                       RETCAP = 15 )
ParDFlow1 = data.frame( B_SOIL = 0.03, C_MAX = 5, B_EVAP = 0.5,  KS = 0.002, KF = 0.2, ADIV = 0.01, CDIV = 0.05,
                        SDIV = 0.01, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.01, TETR = -1, DDFA = 0.08, TMEL = -8.0,
                        RETCAP = 2 )

ParDF2 = data.frame( B_SOIL = 1.6, C_MAX = 500, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                     SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                     RETCAP = 2 )

ParDFup2 = data.frame( B_SOIL = 2, C_MAX = 200, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.3,
                       SDIV = 0.3, CAN_ST = 4., STEM_ST = 4., CSDIV = 0.8, TETR = 0.5, DDFA = 9, TMEL = 0.0,
                       RETCAP = 15 )
ParDFlow2 = data.frame( B_SOIL = 0.03, C_MAX = 5, B_EVAP = 0.5,  KS = 0.002, KF = 0.2, ADIV = 0.01, CDIV = 0.05,
                        SDIV = 0.01, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.01, TETR = -1, DDFA = 0.08, TMEL = -8.0,
                        RETCAP = 2 )
nPars1HRU= ncol(ParDFlow1)


dny=c(30,60,90,120,150,180,210,240,270,300,330,355,364)
p_OBS=dny/365.25

RaKL = 96# odhad Martin Hanel
# ZHU Qa
Qa=10
# mdenni Q
QmKL = c(22, 15, 12, 10, 8.5, 6.5, 6.0, 5.0, 3.5, 3.0, 2.0, 1.0, 0.5)
A=3.28*1000*1000# plocha KL
RmKL = QmKL * (3600*24) / A #CHMU ZHU
Ra =Qa*3600*24*365/A
# simRM=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS)))
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
                    SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )

mae = function(myPar){
  # parrvv=as.numeric(c(ParDFlow1[1,],ParDFlow2[1,]))
  # Par1 = parrvv[1:nPars1HRU]
  # Par2 = as.numeric(parrvv[(nPars1HRU+1):(2*nPars1HRU)])
  # Par1 =ParDF1[1,]
  # Par2 =ParDF2[1,]
  Par1 = myPar[1:nPars1HRU]
  Par2 = as.numeric(myPar[(nPars1HRU+1):(2*nPars1HRU)])
  setParamsToOnedHru(dHRUM_ptr = dhrus,as.numeric(Par1),names(ParDF),0)
  setParamsToOnedHru(dHRUM_ptr  = dhrus,as.numeric(Par2),names(ParDF),1)
  # # for( i in 1:1000){
  calcHBInAlldHrus(dHRUM_ptr = dhrus)
  gatherHBdata(dHRUM_ptr = dhrus)
  # # }
  dta <- getOutput(dHRUM_ptr = dhrus)
  dF <-data.frame(dta$outDta)
  names(dF) <- dta$VarsNams
  simRM=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
  mymae =as.double(sum(abs(simRM - RmKL)))
  # # mymae=NA
  # if(is.na(mymae)) mymae = 9999
  # (return as.double(mymae))
}

library(RcppDE)
itermaxW=10
decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 300,
                itermax = itermaxW, CR = 0.35, F = 0.58, trace = TRUE,
                initialpop = NULL, storepopfrom = itermaxW + 1,
                storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                steptol = itermaxW)

u=DEoptim( lower=as.numeric(c(ParDFlow1[1,],ParDFlow2[1,])),
           upper=as.numeric(c(ParDFup1[1,],ParDFup2[1,])), fn=mae, control = decntr)

u$optim$bestmem

ParBest = as.numeric(u$optim$bestmem)
# ParDF
setParamsToOnedHru(dHRUM_ptr = dhrus,as.numeric(ParBest[1:nPars1HRU]),names(ParDF),0)
setParamsToOnedHru(dHRUM_ptr  = dhrus,as.numeric(ParBest[(nPars1HRU+1):(2*nPars1HRU)]),names(ParDF),1)


calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
# # }
dta <- getOutput(dHRUM_ptr = dhrus)
dF <-data.frame(dta$outDta)

names(dF) <- dta$VarsNams

plot(dF$TOTR, type='l')
plot(dF$BASF, type='l')
plot(dF$DIRR, type='l')
plot(dF$SOIS, type='l')
plot(dF$GROS, type='l')
plot(dF$SNOW, type='l')
plot(dF$MELT, type='l')
simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))


plot(RmKL, simBest)

plot(p_OBS,RmKL, pch=19, main = "Karlův Luh - kalibrace dHRUM - 2 HRU",xlab="P(Qm)", ylab="Qm [mm/den]")
points(p_OBS,simBest,col="red",pch=19)
grid()

ParBest[1:nPars1HRU]
ParBest[(nPars1HRU+1):(2*nPars1HRU)]


require(data.table)
pars=data.table(ParDF1)

ParDF11=ParDF1
ParDF11[1,] = ParBest[1:nPars1HRU]

ParDF11=rbind(ParDF11,ParBest[(nPars1HRU+1):(2*nPars1HRU)])

pars2hru_kl <-cbind(ParDF11,Areas=Areas)

saveRDS(pars2hru_kl,file = "./data/Params/Amalie/KL_2HRU.rds")
