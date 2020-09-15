library(dHRUM)
library(hydroGOF)
library(RcppDE)

# 01	01181000	     WEST BRANCH WESTFIELD RIVER AT HUNTINGTON, MA	  42.23731	 -72.89565	    243.50

dta<-read.table("/home/hubert/Desktop/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/basin_mean_forcing/daymet/01/01181000_lump_cida_forcing_leap.txt",skip=4)
dtaQ<-read.table("/home/hubert/Desktop/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/01/01181000_streamflow_qc.txt")

pr<-dta$V6
temp<-(dta$V9+dta$V10)/2

area <- 243.50 * 1000 * 1000

lat <- 42.23731

Qm <-dtaQ$V5 * 0.0283168466 * 3600 * 24/area*1000


nt=6000
sum(pr[1:nt])
sum(Qm[1:nt])

nHrus <- 1
IdsHrus <- paste0("ID",01181000)
dhrus <- initdHruModel(nHrus,area,IdsHrus)

setPTInputsToAlldHrus(dhrus, Prec = pr[1:nt], Temp = temp[1:nt], inDate = as.Date("1980/01/01"))
calcPetToAllHrus(dHRUM_ptr = dhrus,Latitude = lat,"Hamon")

ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 10 )

ParDFup = data.frame( B_SOIL = 2, C_MAX = 200, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.3,
                      SDIV = 0.3, CAN_ST = 4., STEM_ST = 4., CSDIV = 0.8, TETR = 0.5, DDFA = 10, TMEL = 0.0,
                      RETCAP = 15 )
ParDFlow = data.frame( B_SOIL = 1.3, C_MAX = 5, B_EVAP = 0.5,  KS = 0.002, KF = 0.2, ADIV = 0.01, CDIV = 0.05,
                       SDIV = 0.01, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.01, TETR = -1, DDFA = 0.08, TMEL = -8.0,
                       RETCAP = 2 )

ParBest = ParDF

fitness = function(myPar){
  # myPar =ParDF[1,]
  setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(myPar),names(ParDF))
  # # for( i in 1:1000){
  calcHBInAlldHrus(dHRUM_ptr = dhrus)
  gatherHBdata(dHRUM_ptr = dhrus)
  # # }
  dta <- getOutput(dHRUM_ptr = dhrus)
  dF <-data.frame(dta$outDta)
  names(dF) <- dta$VarsNams
  # print((-1)*as.double(KGE(sim=dF$TOTR[365:nt],obs=Qm[365:nt])))
  sim=dF$TOTR[365:nt]
  obs=Qm[365:nt]
  myfitnes = sum((sim-obs)^2) / sum((obs- mean(obs))^2)
  # # mymae=NA
  # if(is.na(mymae)) mymae = 9999
  # (return as.double(mymae))
}

itermaxW=40

decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 200,
                        itermax = itermaxW, CR = 0.75, F = 0.9, trace = TRUE,
                        initialpop = NULL, storepopfrom = itermaxW + 1,
                        storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                        steptol = itermaxW)

u=DEoptim( lower=as.numeric(ParDFlow[1,]), upper=as.numeric(ParDFup[1,]), fn=fitness, control = decntr)
u$optim$bestmem
ParBest[1,] = as.numeric(u$optim$bestmem)
ParBest

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


plot(Qm[1:nt],type="l")
lines(dF$TOTR,col="red")


KGE(sim=dF$TOTR[365:nt],obs=Qm[365:nt])
NSE(sim=dF$TOTR[365:nt],obs=Qm[365:nt])

plot(Qm[1:nt],dF$TOTR)
