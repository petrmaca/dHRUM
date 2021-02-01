library(dHRUM)
library(hydroGOF)
library(RcppDE)
library(data.table)

# 01	01181000	     WEST BRANCH WESTFIELD RIVER AT HUNTINGTON, MA	  42.23731	 -72.89565	    243.50


pathToCamel <- "/home/hubert/prg/data/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2"
pathToForcing <- "basin_mean_forcing/maurer/"
# pathToForcing <-"/basin_mean_forcing/daymet/"
# pathToForcing <- "basin_mean_forcing/nldas/"
pathToObsQ <- "/usgs_streamflow/"

basinChrs <- fread(paste0(paste0(pathToCamel,"/basin_metadata/basin_physical_characteristics.txt")))
gaugeChars <- fread(paste0(paste0(pathToCamel,"/basin_metadata/gauge_informationPM.txt")))

gaugeChars[523:524,]

i <- 1
nt <- 2000
nseVec <- c()
kgeVec <- c()
for(i in 1:nrow(gaugeChars)){
  ifelse((gaugeChars$GAGE_ID[i]<10000000), gChrs <-paste0(0,gaugeChars$GAGE_ID[i]), gChrs <- paste0(gaugeChars$GAGE_ID[i]))
  ifelse((gaugeChars$HUC_02[i]<10), HUc <-paste0(0,gaugeChars$HUC_02[i],"/"), HUc <- paste0(gaugeChars$HUC_02[i],"/"))
  dtaFC <- read.table(paste0(pathToCamel,pathToForcing,HUc,gChrs,"_lump_cida_forcing_leap.txt"), skip = 4)
  # dtaFC <- read.table(paste0(pathToCamel,pathToForcing,HUc,gChrs,"_lump_maurer_forcing_leap.txt"), skip = 4)
  # dtaFC <- read.table(paste0(pathToCamel,pathToForcing,HUc,gChrs,"_lump_nldas_forcing_leap.txt"), skip = 4)
  dtaQ <- read.table(paste0(pathToCamel,pathToObsQ,HUc,gChrs,"_streamflow_qc.txt"))

  pr<-dtaFC$V6
  temp<-(dtaFC$V9+dtaFC$V10)/2

  area <- gaugeChars$`DRAINAGE_AREA_(KM^2)`[i] * 1000 * 1000
  lat <- gaugeChars$LAT[i]

  Qm <-dtaQ$V5 * 0.0283168466 * 3600 * 24/area*1000

  print(paste("The basin number: ",i))
  print(sum(pr[1:nt]))
  print(sum(Qm[1:nt]))

  nHrus <- 1
  IdsHrus <- paste0("ID_",0,gaugeChars$GAGE_ID[i])
  dhrus <- initdHruModel(nHrus,area,IdsHrus)

  setPTInputsToAlldHrus(dhrus, Prec = pr[1:nt], Temp = temp[1:nt], inDate = as.Date("1980/01/01"))
  calcPetToAllHrus(dHRUM_ptr = dhrus,Latitude = lat,"Hamon")

  ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
                    SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
                    RETCAP = 10 )
  ParDFup = data.frame( B_SOIL = 2, C_MAX = 800, B_EVAP = 2,  KS = 0.4, KF = 0.7, ADIV = 0.9, CDIV = 0.3,
                      SDIV = 0.3, CAN_ST = 4., STEM_ST = 4., CSDIV = 0.8, TETR = 0.5, DDFA = 10, TMEL = 0.0,
                      RETCAP = 25 )
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

  itermaxW=100
  decntr<-DEoptim.control(VTR = 0, strategy = 2, bs = FALSE, NP = 300,
                        itermax = itermaxW, CR = 0.75, F = 0.9, trace = FALSE,
                        initialpop = NULL, storepopfrom = itermaxW + 1,
                        storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
                        steptol = itermaxW)
  u=DEoptim( lower=as.numeric(ParDFlow[1,]), upper=as.numeric(ParDFup[1,]), fn=fitness, control = decntr)
  ParBest[1,] = as.numeric(u$optim$bestmem)

  setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParBest[1,]),names(ParDF))

  calcHBInAlldHrus(dHRUM_ptr = dhrus)
  gatherHBdata(dHRUM_ptr = dhrus)

  dta <- getOutput(dHRUM_ptr = dhrus)
  dF <-data.frame(dta$outDta)
  names(dF) <- dta$VarsNams

  kgeVec[i] <- KGE(sim=dF$TOTR[365:nt],obs=Qm[365:nt])
  nseVec[i] <- NSE(sim=dF$TOTR[365:nt],obs=Qm[365:nt])
  plot(Qm[1:nt],dF$TOTR, main=paste(gaugeChars$GAGE_ID[i]," KGE=", format(kgeVec[i], digits = 2)," NSE=", format(nseVec[i], digits = 2)))
  plot(Qm[1:nt],type="l", main=paste(gaugeChars$GAGE_ID[i]," KGE=", format(kgeVec[i], digits = 2)," NSE=", format(nseVec[i], digits = 2))
       ,ylab ="Q [mm]", xlab="Time [day]")
  lines(dF$TOTR,col="red")
  grid()
}

boxplot(data.frame(KGE =kgeVec,NSE=nseVec))
KGE_NSE=data.frame(KGE =kgeVec,NSE=nseVec)
saveRDS(file="/home/hubert/prg/data/outputs/NT_2000/camel_kge_nse_basic_2000cal_Daymet.rds",KGE_NSE)
# saveRDS(file="/home/hubert/prg/data/outputs/NT_2000/camel_kge_nse_basic_2000cal_Maurer.rds",KGE_NSE)
# saveRDS(file="/home/hubert/prg/data/outputs/NT_2000/camel_kge_nse_basic_2000cal_nldas.rds",KGE_NSE)

# dtKGdaymet <- readRDS("/home/hubert/prg/data/outputs/NT_2000/camel_kge_nse_basic_2000cal_Daymet.rds")
# dtKGmauer <- readRDS("/home/hubert/prg/data/outputs/NT_2000/camel_kge_nse_basic_2000cal_Maurer.rds")
# dtKGnldas <- readRDS("/home/hubert/prg/data/outputs/NT_2000/camel_kge_nse_basic_2000cal_nldas.rds")

KGE = data.frame(KGEdaymet = dtKGdaymet$KGE, KGEmaurer = dtKGmauer$KGE, KGEnldas = dtKGnldas$KGE)
boxplot(KGE)
summary(KGE)

NSE = data.frame(NSEdaymet = dtKGdaymet$NSE, NSEmaurer = dtKGmauer$NSE, NSEnldas = dtKGnldas$NSE)
boxplot(NSE)
summary(NSE)

#
# plot(dF$TOTR, type='l')
# plot(dF$BASF, type='l')
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
#
#
# plot(Qm[1:nt],type="l")
# lines(dF$TOTR,col="red")



# plot(Qm[1:nt],dF$TOTR)
