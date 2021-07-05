library(dHRUM)
library(hydroGOF)
library(RcppDE)
library(data.table)
library(rgdal)
library(readxl)

dtaFC = readRDS("~/CULS_FES/dHRUM/TACR_Dyje/meteo_hamr.rds")
dtaFC <- dtaFC[order(DBC),]

#dtaFC[21848]
#dtaFC[43695]
dtaFC[87389]

DBC = as.numeric(dtaFC[,unique(DBC),])
length(DBC)
DBC

pathToObsQ <- "~/CULS_FES/dHRUM/TACR_Dyje/Hydrologicka_data/data/"
pathToObsArea <- "~/CULS_FES/dHRUM/TACR_Dyje/povodi/"
#dtaFC <- dtaFCc[87389:length(dtaFCc)]

files <- list.files(pathToObsQ)

for(i in 1:length(files)) {
  files[i] = substr(files[i],1,4)
}

# files = as.numeric(files)
# nov <- setdiff(DBC,files)
# nov[1]
# DBC[- c(4276)]
#
# DBC
# for(i in 1:length(files)){
#   if(files[i] %in% nov)
# }


i <- 2
nt <- 2000
nseVec <- c()
kgeVec <- c()
for(i in 2:length(DBC)) {
  basinId <- paste0(as.character(DBC[i]))

  dtaQ <- fread(paste0(paste0(pathToObsQ,basinId,".txt")))

  gisArea <- readOGR(
    dsn = paste0(pathToObsArea),
    layer = basinId,
    verbose = FALSE
  )

  area <- gisArea$AREA

  xlsxLat <- read_excel("~/CULS_FES/dHRUM/TACR_Dyje/Hydrologicka_data/seznam_objektu_pvv_pzv.xlsx", sheet="Dyje - PVV")

  i <- 1
  for(i in 1:nrow(xlsxLat)) {
    num_id <- paste0(basinId, "00")
    num_id <- as.numeric(num_id)
    if (xlsxLat$DBC[i] == num_id) {
      lat <- xlsxLat$Y_WGS[i]
    }
  }

  dtaQ$Q <- as.numeric(sub(",", ".", sub(".", "", dtaQ$Q, fixed=TRUE), fixed=TRUE))
  Qm <- dtaQ$Q * 3600 * 24 / area * 1000

  pr <- dtaFC$PR
  temp <- dtaFC$TAS

  pr <- pr[25500:27499]
  temp <- temp[25500:27499]

  print(sum(Qm[1:nt]))
  print(sum(pr[1:nt]))

  dtaFC$DTM <- as.Date(dtaFC$DTM,
                         format = "%Y/%m/%d")

  nHrus <- 1
  dhrus <- initdHruModel(nHrus,area,basinId)

  setPTInputsToAlldHrus(dhrus, Prec = pr[1:nt], Temp = temp[1:nt], inDate = as.Date("1971/01/01"))
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
                          itermax = itermaxW, CR = 0.75, F = 0.9, trace = TRUE,
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
  plot(Qm[1:nt],dF$TOTR, main=paste(basinId," KGE=", format(kgeVec[i], digits = 2)," NSE=", format(nseVec[i], digits = 2)))
  plot(Qm[1:nt],type="l", main=paste(basinId," KGE=", format(kgeVec[i], digits = 2)," NSE=", format(nseVec[i], digits = 2))
       ,ylab ="Q [mm]", xlab="Time [day]")
  lines(dF$TOTR,col="red")
}
