library(dHRUM)
library(lubridate)
library(data.table)
nHrus <- 1
Areas <- runif(nHrus,min = 1,max  = 100)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)

filname2 = "../tests/indata/KL_1960_01_01.txt"
TPdta = read.table(filname2)

prec=TPdta$V1
temp=TPdta$V2
setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, inDate = as.Date("1960/01/01"))
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Hamon")

ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 10 )

ParDFup = data.frame( B_SOIL = 2, C_MAX = 200, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.3,
                      SDIV = 0.3, CAN_ST = 4., STEM_ST = 4., CSDIV = 0.8, TETR = 0.5, DDFA = 5, TMEL = 0.0,
                      RETCAP = 15 )
ParDFlow = data.frame( B_SOIL = 0.03, C_MAX = 5, B_EVAP = 0.5,  KS = 0.002, KF = 0.2, ADIV = 0.01, CDIV = 0.05,
                       SDIV = 0.01, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.01, TETR = -1, DDFA = 0.08, TMEL = -8.0,
                       RETCAP = 2 )

ParBest = ParDF
  # set_Params_TodHru(dHRU_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF),TRUE,0)
  # # for( i in 1:1000){
  # calc_HBInAlldHrus(dHRU_ptr = dhrus)
  # gatherHBdata(dHRU_ptr = dhrus)
  # # }
  # dta <- get_Output(dHRU_ptr = dhrus)
  # dF <-data.frame(dta$outDta)
  # names(dF) <- dta$VarsNams
  # dta

dny=c(30,60,90,120,150,180,210,240,270,300,330,355,364)
p_OBS=dny/365.25

# BP 96.68961 (1.736-2.430446)
RaKL = 96# odhad Martin Hanel
# ZHU Qa
Qa=10
# mdenni Q
QmKL = c(22, 15, 12, 10, 8.5, 6.5, 6.0, 5.0, 3.5, 3.0, 2.0, 1.0, 0.5)
A=3.28*1000*1000# plocha KL
RmKL = QmKL * (3600*24) / A #CHMU ZHU
Ra =Qa*3600*24*365/A
# simRM=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS)))

mae = function(myPar){
 # myPar =ParDF[1,]
  setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(myPar),names(ParDF))
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
itermaxW=15

decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 200,
                itermax = itermaxW, CR = 0.75, F = 0.9, trace = TRUE,
                initialpop = NULL, storepopfrom = itermaxW + 1,
                storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                steptol = itermaxW)
n_ens=25
parsKLmatrix=matrix(0,nrow=n_ens, ncol=ncol(ParBest))
for(i in 1:n_ens){
  
u=DEoptim( lower=as.numeric(ParDFlow[1,]), upper=as.numeric(ParDFup[1,]), fn=mae, control = decntr)

u$optim$bestmem

# ParBest[1,] = as.numeric(u$optim$bestmem)
# ParBest
# u$optim$bestmem
ParBest[1,] = as.numeric(u$optim$bestmem)
parsKLmatrix[i,] = as.numeric(u$optim$bestmem)
}

KL_df = data.frame(parsKLmatrix)
names(KL_df) = names(ParBest)
KL_df=cbind(KL_df,ID=rep("KL",times=n_ens))

A=3.28*1000*1000# plocha KL
Par_dHRUm_KL_lumped= list(
  pars =KL_df,
  areaM2 = A,
  Pet_type = "Hamon",
  latitude =50.1
)
save(Par_dHRUm_KL_lumped,file="par_dHRUM_lumped_PET_HAMON_50_1_KL.rda")




# ParDF
setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParBest[1,]),names(ParDF))
# for( i in 1:1000){
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
# }
dta <- getOutput(dHRUM_ptr = dhrus)
dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams

plot(dF$TOTR, type='l')
plot(dF$BASF, type='l')
plot(dF$DIRR, type='l')
plot(dF$SOIS, type='l')
plot(dF$SURS,type= 'l')
plot(dF$GROS, type='l')
plot(dF$PREF,type= 'l')
plot(dF$TROF,type= 'l')
plot(dF$SNOW, type='l')

plot(dF$CANS, type='l')
plot(dF$STES, type='l')

plot(dF$EVAC, type='l')
plot(dF$EVAS, type='l')
plot(dF$EVBS, type='l')
plot(dF$AET, type='l')

which.min(dF$EVBS)

aet=dF$AET+dF$EVBS+dF$EVAC+dF$EVAS
aet1=dF$AET+dF$EVBS
plot(aet)

plot((dF$AET+dF$EVBS+dF$EVAC+dF$EVAS)[1:700],type='l')
plot((dF$SOIS)[1:1700],type='l')

dFF= data.table(
  DTM = as.Date(paste0(dF$YEAR,"-",dF$MONTH,"-",dF$DAY)),
  GW = dF$GROS,
  SO = dF$SOIS,
  SR = dF$SURS,
  P=dF$PREC,
  Sn=dF$SNOW,
  T=dF$TEMP,
  PET=dF$PET,
  AET=dF$AET+dF$EVBS+dF$EVAC+dF$EVAS,
  MEL =dF$MELT,
  R =dF$TOTR
)
dFF[,Month:=month(DTM)]
dFF[,Year:=year(DTM)]

dFF[DTM>as.Date("2006-01-01"),.(GW=mean(GW),SO=mean(SO),SR=mean(SR), SOSR=mean(SO+SR),P=sum(P),Sn=mean(Sn),T=mean(T),PET=sum(PET),AET=sum(AET),MELt=sum(MEL),R=sum(R)),by=.(Month)]
dFF[DTM>as.Date("1999-01-01"),.(GW=mean(GW),SO=mean(SO),SR=mean(SR), SOSR=mean(SO+SR),P=sum(P),Sn=mean(Sn),T=mean(T),PET=sum(PET),AET=sum(AET),MELt=sum(MEL),R=sum(R)),by=.(Year)]
simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))

simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
  plot(RmBP, simBest)

plot(p_OBS,RmBP, pch=19)
points(p_OBS,simBest,col="red",pch=19)

saveRDS(file="../dHRUM/data/basin_params/KL/dHRUM_lumped_KL_01.rds",ParBest)

# readRDS("../tests/param0dKL.rds")