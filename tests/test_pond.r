# library(profvis)
# profvis(
# {
library(dHRUM)
library(data.table)

#PRO TVORBU GRAFU
library(dygraphs)
library(xts)
library(htmlwidgets)

numdata =365
probwet =0.9
meanifwet = 8

prec= rbinom(numdata,1,probwet)*rexp(numdata,1/meanifwet)
temp=rnorm(numdata,20,3)
lai = abs(rnorm(numdata,3,2))

par_values=c()
par_names=c()

#nHrus <- 15000
nHrus <- 1
Areas <- runif(nHrus,min = 38780000,max  = 38780050) #range of areas for test pond situation


IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes=rep("DDF",times= length(Areas)),hruIds=IdsHrus)

setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus,
                              InstStLai=rep(TRUE,times= length(Areas)),
                              smaxlaiTypes = rep("VonHoyningenHuene",times= length(Areas)))

setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
setFastSlowDivModeltypeToAlldHrus(dHRUM_ptr = dhrus, adivsModelTypes = rep("CnstAdiv",times= length(Areas)), hruIds = IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))






################################ test pond - start #############################################################################
#setPondToAlldHrus(dHRUM_ptr = dhrus,PondTypes=rep("Pond",times= length(Areas)),hruIds=IdsHrus)
# pondDF1 and pondDF2 are pond variables and specifications
pondDF1 = data.frame( PondArea = 40500, PonsMax= 45000, MRF= 0.039, Coflw=0.3)
pondDF2 = data.frame( Pond_ET = "ETpond1", Pond_inSOIS= "noPondSOISPerc", Pond_inGW = "noPondGWPerc",
                      Pond_outSOIS= "noPondSOISPerc", Pond_outGW= "noPondGWPerc",Pond_outReg="PondRouT3" )
# pond implementation
setPondToOnedHru(dHRUM_ptr = dhrus,0,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
################################ test pond - end #############################################################################


calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
outDta <- dHRUMrun(dHRUM_ptr = dhrus)

outDF <- data.frame(outDta$outDta)
names(outDF) <-c(outDta$VarsNams)
outDF = as.data.table(outDF)

################################ pond - visualisation #############################################################################
plot(outDF$PONS, type ="l",ylim=c(0,60000),main = "Response of water volume in the dam to precipitation",xlab="day")
lines(60000-(outDF$PERC*300), type ="l",col="RED")
legend("bottomright", legend = c("PERC*300 ", "PONS"),
       text.width = strwidth("1,000,000"),
       lty = c(1,1), xjust = 1, yjust = 1,col=c("red","black"),
       title = "Line Types")



DTM=as.Date(paste(outDF$YEAR, outDF$MONTH, outDF$DAY, sep = "-"))
no_pond <- xts(x = outDF$TOTR, order.by = DTM)
# pond <- xts(x = outDF_pond$TOTR, order.by = DTM)

# don <- cbind(no_pond,pond)
don <- cbind(no_pond)
mmain=paste0("Graf simulace odtoku z povodi ")
p <- dygraph(don,main = mmain) %>%
  dyAxis("y", label = "Průměrný denní odtok v mm") %>%
  dyOptions(labelsUTC = TRUE, fillGraph=TRUE, fillAlpha=0.1, drawGrid = FALSE, colors = c("BLACK","RED","BLUE","GREEN")) %>%
  dyRangeSelector() %>%
  dyCrosshair(direction = "vertical") %>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE)  %>%
  dyRoller(rollPeriod = 1)
p

