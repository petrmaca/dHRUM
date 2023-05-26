library(dHRUM)
library(RcppDE)
library(hydroGOF)
library(data.table)

# numdata =14671
# probwet =0.3
# meanifwet = 0.5
# 
# namsP= dir("./data/meteo156stanic/P/")
# 
# data.table(readRDS("./data/meteo156stanic/P/P_1989.rds"))
# data.table(readRDS("./data/meteo156stanic/T/T_1989.rds"))
# data.table(readRDS("./data/meteo156stanic/E/EVAP_1969.rds"))
# 
# 
# length(namsP)
# prec=list()
# for(i in 1:2){
# prec[[i]]=data.table(readRDS(paste0("./data/meteo156stanic/P/",namsP[20])))
# }
# prec=rbindlist(prec)
# prec[ID %in% 905,]
# # 
# prec= rbinom(numdata,1,probwet)*rexp(numdata,1/meanifwet)
# temp=20*sin((1:14671)/45+35)+rnorm(numdata,0,2)
# plot(temp,type="l")

# dtaMPX <- read.table("./data/mopex/12027500.txt")


dtaMPX <- read.table("/home/hubert/Desktop/terka/data/mopex/12027500g.csv",header=TRUE)


# nrow(dtaMPX[dtaMPX$V1==1954,])
names(dtaMPX)=c("Year","Month","Day","Pre","Pet","R","Tmax","Tmin")
plot(dtaMPX$R[19996:19997])
dtaMPX=as.data.table(dtaMPX)
dtaMPX[,T:=(Tmax+Tmin)/2]


# dF
# dF$BASF

# plot(20*sin((1:365)/45+35), type="l")


# kalibrace
# ParDF = data.frame( B_SOIL = 1, C_MAX = 10, CMIN=5, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
#                     SDIV = 0.1, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = -1.0,
#                     RETCAP = 2, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
#                     KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3)
# 
# ParDFup = data.frame( B_SOIL = 2, C_MAX = 200, CMIN=5, B_EVAP = 1,  KS = 0.1, KF = 0.3, ADIV = 0.8, CDIV = 0.3,
#                       SDIV = 0.1, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = -1.0,
#                       RETCAP = 2, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
#                       KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3)
# 
# ParDFlow = data.frame(B_SOIL = 1, C_MAX = 1, CMIN=5, B_EVAP = 1,  KS = 0.001, KF = 0.003, ADIV = 0.8, CDIV = 0.3,
#                       SDIV = 0.1, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = -1.0,
#                       RETCAP = 2, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
#                       KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3)

nHrus <- 1
area= 1
nt =5000
Qm=dtaMPX$R[1:nt]
temp=dtaMPX$T
pr=dtaMPX$Pre
lat=45

plot(Qm)

IdsHrus <- paste0("ID_",0,"leafriver")
dhrus <- initdHruModel(nHrus,area,IdsHrus)
gwStorType <- c("LINL_RES")
swStorType <- c("PDM2")
setGWtypeToAlldHrus(dHRUM_ptr = dhrus, gwStorType, IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, swStorType, IdsHrus)
setPTInputsToAlldHrus(dhrus, Prec = pr[1:nt], Temp = temp[1:nt], inDate = as.Date("1980/01/01"))
calcPetToAllHrus(dHRUM_ptr = dhrus,Latitude = lat,"HAMON")

ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 10, D_BYPASS = 0.8, THR = 10, KS2 = 0.1, ALPHA = 0.5, FOREST_FRACT = 0.3, FC = 10,
                    KF_NONLIN = 10, KF2 = 0.01, C = 10, INFR_MAX = 10, RF = 0.5, WP = 0.3,CMIN =20,L=0.1)
ParDFup = data.frame( B_SOIL = 2, C_MAX = 800, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.3,
                      SDIV = 0.3, CAN_ST = 4., STEM_ST = 4., CSDIV = 0.8, TETR = 0.5, DDFA = 10, TMEL = 0.0,
                      RETCAP = 25, D_BYPASS = 1, THR = 100, KS2 = 0.4, ALPHA = 1, FOREST_FRACT = 0.3, FC = 100,
                      KF_NONLIN = 100, KF2 = 0.1, C = 100, INFR_MAX = 100, RF = 1, WP = 14,CMIN =20,L=2.0)
ParDFlow = data.frame( B_SOIL = 1.3, C_MAX = 50, B_EVAP = 0.5,  KS = 0.002, KF = 0.2, ADIV = 0.01, CDIV = 0.05,
                       SDIV = 0.01, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.01, TETR = -1, DDFA = 0.08, TMEL = -8.0,
                       RETCAP = 2, D_BYPASS = 0.1, THR = 1, KS2 = 0.002, ALPHA = 0.1, FOREST_FRACT = 0.01, FC = 1,
                       KF_NONLIN = 1, KF2 = 0.001, C = 1.0, INFR_MAX = 1.0, RF = 0.01, WP = 0.30,CMIN =0,L=-2.0)

ParBest = ParDF



fitness = function(myPar){
  # myPar =ParDF[1,]
  
  setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(myPar),names(ParDF))
  # # for( i in 1:1000){
  outDta <- dHRUMrun(dHRUM_ptr = dhrus)
  outDF <- data.frame(outDta$outDta)
  names(outDF) <- outDta$VarsNams
  # # }
  dF <- outDF
  # print((-1)*as.double(KGE(sim=dF$TOTR[365:nt],obs=Qm[365:nt])))
  sim=log(dF$TOTR[365:nt])
  obs=log(dtaMPX$R[365:nt])
  
  mae = mae(sim = sim, obs = obs)
  
  if (is.na(mae)) {
    mae = 99999999
  }
  
  mae
  #myfitnes = sum((sim-obs)^2) / sum((obs- mean(obs))^2)
  # # mymae=NA
  # if(is.na(mymae)) mymae = 9999
  # (return as.double(mymae))
}
# n_ens=1
# parsKLmatrix=matrix(0,nrow=n_ens, ncol=ncol(ParBest))
# for(i in 1:n_ens){

itermaxW=10
decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 300,
                        itermax = itermaxW, CR = 0.75, F = 0.9, trace = TRUE,
                        initialpop = NULL, storepopfrom = itermaxW + 1,
                        storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                        steptol = itermaxW)
u=DEoptim( lower=as.numeric(ParDFlow[1,]), upper=as.numeric(ParDFup[1,]), fn=fitness, control = decntr)

u$optim$bestmem

ParBest[1,] = as.numeric(u$optim$bestmem)
ParBest
u$optim$bestmem

setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParBest[1,]),names(ParDF))
#run in a for loop and update state variables odnosno reset za sekoja presmetka
#preserve mass balance
outDta <- dHRUMrun(dHRUM_ptr = dhrus)
outDF <- data.frame(outDta$outDta)
names(outDF) <- outDta$VarsNams
dF <- outDF
dF
plot(outDF$SOIS)
plot(outDF$PERC)
plot(outDF$GROS)
plot(Qm[1:nt],dF$TOTR[1:nt], main=paste(" KGE=", format(KGE(sim=dF$TOTR[365:nt],obs=dtaMPX$R[365:nt]), digits = 2)," NSE=", format(NSE(sim=dF$TOTR[365:nt],obs=dtaMPX$R[365:nt]), digits = 2)))
plot(Qm[1:nt],type="l", main=paste(" KGE=", format(KGE(sim=dF$TOTR[365:nt],obs=dtaMPX$R[365:nt]), digits = 2)," NSE=", format(NSE(sim=dF$TOTR[365:nt],obs=dtaMPX$R[365:nt]), digits = 2))
     ,ylab ="Q [mm]", xlab="Time [day]")
lines(dF$TOTR[1:nt],col="red")
grid()
