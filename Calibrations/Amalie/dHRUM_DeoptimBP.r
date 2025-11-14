# sessionInfo()
library(dHRUM)
library(data.table)
nHrus <- 1
Areas <- 4.7*1000*1000
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)


setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_RES",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus)

filname2 = "../dHRUM/Calibrations/Amalie/indata/BP_1960_01_01.txt"

setPTInputsToAlldHrusFromFile(dHRUM_ptr = dhrus, filname2)
calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")

ParDF = data.frame( B_SOIL = 1.6, C_MAX = 35, B_EVAP = 2.5,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.2,
                    SDIV = 0.1, CAN_ST = 2, STEM_ST = 1, CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 10, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
                    KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3,CMIN =25,L=0.1, B_EXP = 0.3, KFR = 0.03)


# ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
                    # SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(ParDF[1,]),ParsNames =names(ParDF))
myPars=getCurdHRUpars(dHRUM_ptr = dhrus,0)

ups =names(ParDF)

myPars

ups1 =  myPars$Cur_names

setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(myPars[[2]]),ParsNames =ups1)
myPars=getCurdHRUpars(dHRUM_ptr = dhrus,0)

# ups =names(ParDF)

myPars

# getCurdHRUpars(dHRUM_ptr = dhrus,0)
# ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
#                     SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
#
# setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(ParDF[1,]),ParsNames =names(ParDF))

# ParDF = data.frame( B_SOIL = 1, C_MAX = 100, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.01, CDIV = 0.3,
# SDIV = 0.3, CAN_ST = .5, STEM_ST = .2, CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
# RETCAP = 2,CMIN =10)
#
# ParDFup = data.frame( B_SOIL = 1, C_MAX = 200, B_EVAP = 2,  KS = 0.8, KF = 0.08, ADIV = 0.99, CDIV = 0.3,
#                       SDIV = 0.3, CAN_ST = 0.5, STEM_ST = .1, CSDIV = 0.8, TETR = 0.5, DDFA = 5, TMEL = 0.0,
#                       RETCAP = 6 ,CMIN =5)
#
# ParDFlow = data.frame( B_SOIL = 1, C_MAX = 21, B_EVAP = 1,  KS = 0.001, KF = 0.005, ADIV = 0.71, CDIV = 0.05,
#                        SDIV = 0.01, CAN_ST = 0.02, STEM_ST = 0.01, CSDIV = 0.01, TETR = -1, DDFA = 0.08, TMEL = -1.0,
#                        RETCAP = 2 ,CMIN =1)
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
RaBP = 96# odhad Martin Hanel
QmBP = c(26, 18, 14, 12, 10, 8.0, 7.0, 6.0, 4.5, 3.5, 2.5, 1.0, 0.5)#l/s in 1 day
A=4.7*1000*1000# plocha BP
RmBP = QmBP * (3600*24) / A #CHMU ZHU mm/day
# simRM=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS)))

sse = function(myPar){
   # myPar =ParDF[1,]
  setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(myPars[[2]]),names(myPars[[1]]))
  # # for( i in 1:1000){
  calcHBInAlldHrus(dHRUM_ptr = dhrus)
  gatherHBdata(dHRUM_ptr = dhrus)
  # # }
  dta <- getOutput(dHRUM_ptr = dhrus)
  dF <-data.frame(dta$outDta)
  names(dF) <- dta$VarsNams
  simRM=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
  # mymae =as.double(sum(abs(simRM - RmBP)))
  mymae =as.double(sum((simRM - RmBP)^2))
  # # mymae=NA
  # if(is.na(mymae)) mymae = 9999
  # (return as.double(mymae))
}

library(RcppDE)

itermaxW=6

decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 200,
                itermax = itermaxW, CR = 0.85, F = 0.95, trace = TRUE,
                initialpop = NULL, storepopfrom = itermaxW + 1,
                storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                steptol = itermaxW)

n_ens=1
parsBPmatrix=matrix(0,nrow=n_ens, ncol=ncol(ParBest))
for(i in 1:n_ens){
  u=DEoptim( lower=as.numeric(myPars[[4]]), upper=as.numeric(myPars[[3]]), fn=sse, control = decntr)
  u$optim$bestmem
  ParBest[1,] = as.numeric(u$optim$bestmem)
  parsBPmatrix[i,] = as.numeric(u$optim$bestmem)
}

BP_df = data.frame(parsBPmatrix)
names(BP_df) = names(ParBest)
BP_df=cbind(BP_df,ID=rep("BP",times=n_ens))
# save(BP_df,file="par_dHRUM_lumped_PET_HAMON_50_1_BP.rda")
# save(BP_df,file="./Calibrations/Amalie/outdata/par_dHRUM_lumped_PET_HAMON_50_1_BP.rda")

# load("./Calibrations/Amalie/outdata/par_dHRUM_lumped_PET_HAMON_50_1_BP.rda")
BP_df


A=4.7*1000*1000## plocha BP
Par_dHRUm_BP_lumped= list(
  pars =BP_df,
  areaM2 = A,
  Pet_type = "HAMON",
  latitude =50.1
)
pars=rbind(Par_dHRUm_BP_lumped$pars)

pars=cbind(pars,CMIN=1)

Par_dHRUm_Amalie_lumped = list(
  pars =pars,
  areaM2 = data.frame(BP= 4700000),
  Pet_type = "HAMON",
  latitude =50.1
)

# save(Par_dHRUm_Amalie_lumped ,file = "./Calibrations/Amalie/outdata/Amalie_lumped_dHRUM.rda")
Par_dHRUm_Amalie_lumped$areaM2$BP
BP_df=cbind(ID=rep("BP",times=n_ens))
BP_df=cbind(ID=rep("CMIN",times=n_ens))


ParBest[1,] = as.numeric(u$optim$bestmem)
# ParBest[1,] =BP_df[1,1:15]

setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParBest),names(ParBest))

  # for( i in 1:1000){
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
# }
dta <- getOutput(dHRUM_ptr = dhrus)
dF <-data.frame(dta$outDta)
names(dF) <- dta$VarsNams


dF=as.data.table(dF)
dF[,DM:=paste(DAY,".",MONTH)]

plot(dF$TOTR)

mS=dF[,.(SOIS =mean(SOIS),PR=mean(PREC),SUR = mean(SURS), RET = mean(SURS+SOIS),Ro=mean(TOTR), AEt = mean(AET+EVBS+EVAC+EVAS), PEt =mean(PET)),by=.(DM)]
plot(mS$SOIS,type='l',ylim=c(0,max(mS$SOIS)))
lines(mS$SUR,col="blue")
lines(mS$RET,col="green")
plot(mS$PR)

points(mS$PEt, col="red")
points(mS$AEt)


plot(mS$PEt, col="red")
points(mS$AEt)

plot(mS$Ro)
plot(mS$SOIS)


AAe=(dF$AET+dF$EVBS+dF$EVAC+dF$EVAS)
o=1:1604

plot(dF$SOIS[o])

plot(dF$SOIS)

max(dF$SOIS)
plot(dF$SOIS[o],col='blue',ylim=c(0,max(c(AAe,dF$PET,dF$SOIS))))
lines(dF$EVBS[o],col='darkolivegreen4')
lines(AAe[o],col="red")
lines(dF$PET[o],col='black')


plot(dF$PET[o],type='l')
lines(dF$EVBS[o],col='darkolivegreen4')
lines(AAe[o],col="red")

simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
plot(p_OBS,RmBP, pch=19)
points(p_OBS,simBest,col="red",pch=19)


#
#
#
#
# plot(dF$TOTR, type='l')
# plot(dF$BASF, type='l')
# sum(dF$BASF)/length(dF$BASF)
# plot(dF$DIRR, type='l')
# plot(dF$SOIS, type='l')
# plot(dF$SURS,type= 'l')
# plot(dF$GROS, type='l')
# plot(dF$PREF,type= 'l')
# plot(dF$TROF,type= 'l')
# plot(dF$SNOW, type='l')
#
# plot(dF$CANS, type='l')
# plot(dF$STES, type='l')
#
# plot(dF$EVAC, type='l')
# plot(dF$EVAS, type='l')
# plot(dF$EVBS, type='l')
# plot(dF$AET, type='l')
# plot(dF$MELT, type='l')
#
# which.min(dF$EVBS)
#
# aet=dF$AET+dF$EVBS+dF$EVAC+dF$EVAS
# aet1=dF$AET+dF$EVBS
# plot(aet)
#
#
# # plot((dF$SOIS)[1:1700],type='l')
#
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
  R =dF$TOTR,
  BF = dF$BASF,
  DR = dF$DIRR,
  PR = dF$PERC,
  EF = dF$PREF)
dFF[,Month:=month(DTM)]
dFF[,Year:=year(DTM)]

dFF[DTM>as.Date("2006-01-01"),.(GW=mean(GW),SO=mean(SO),SR=mean(SR), SOSR=mean(SO+SR),P=sum(P),Sn=sum(Sn),T=mean(T),PET=sum(PET),AET=sum(AET),MELt=sum(MEL),R=sum(R), BF=sum(BF), DR=sum(DR), PERC= sum(PR),EF = sum(EF)),by=.(Month)]
HBy =dFF[DTM>as.Date("1961-01-01"),.(GW=mean(GW),SO=mean(SO),SR=mean(SR), SOSR=mean(SO+SR),P=sum(P),Sn=mean(Sn),T=mean(T),PET=sum(PET),AET=sum(AET),MELt=sum(MEL),R=sum(R), BF=sum(BF), DR=sum(DR), PERC= sum(PR)),by=.(Year)]
dFF[DTM>as.Date("1961-01-01"),.(GW=mean(GW),SO=mean(SO),SR=mean(SR), SOSR=mean(SO+SR),P=sum(P),Sn=mean(Sn),T=mean(T),PET=sum(PET),AET=sum(AET),MELt=sum(MEL),R=sum(R)),]

mean(HBy$R)
mean(HBy$PET)
mean(HBy$AET)
# simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
#
# plot((dF$AET+dF$EVBS+dF$EVAC+dF$EVAS)[1:700],type='l')
# plot((dF$SOIS)[1:1700],type='l')
# simBest=as.numeric(quantile(dF$TOTR,probs=(1-p_OBS), na.rm = TRUE))
# plot(RmBP, simBest)
#
# plot(RmBP, simBest)
#
#
#
# plot(p_OBS,RmBP, pch=19)
# points(p_OBS,simBest,col="red",pch=19)
#
# (RmBP-simBest)/RmBP*100
#
# # saveRDS(file="./data/basin_params/BP/dHRUM_lumped_BP_03.rds",ParBest)
#
