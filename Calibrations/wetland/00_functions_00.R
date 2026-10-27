
vyber_model<-function(sel_pov,hlavni_slozka){
  #priklad pouziti funkce
  #vyber_model(sel_pov="1-13-01-0750-0-00-00",hlavni_slozka = "data_od_martina")

#pro rybnik potrebuju nakalibrovanou strukturu a parametry.
#Ty se berou od martina - nejake vypocty, ne vcechny kombinace, ne vsechna povodi, je to nekde na Horce. je to vypocet s behavioralnima metrikama.
#Z martinovych datasetu vytahnu modely pro jedno povodi, ty dle KGE seradima  zkousim vzit nejlepsi

library(data.table)


seznam_souboru <- list.files(
  path = hlavni_slozka,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

dta_1=as.data.table(readRDS(seznam_souboru[1]))

subset_1 <- dta_1[ID == sel_pov]

KGE_max_1 <- subset_1[which.max(OF_TOT)]


return (KGE_max_1)

}


graf_porovnani_TOTR<-function(outDF_ref,outDF_sim_measure,CatchmentArea_m2,uloz_graf){
  DTM=as.Date(paste(outDF_sim_measure$YEAR, outDF_sim_measure$MONTH, outDF_sim_measure$DAY, sep = "-"))

  ref <- xts(x = outDF_ref$TOTR/1000*CatchmentArea_m2, order.by = DTM)
  measure <- xts(x = outDF_sim_measure$TOTR/1000*CatchmentArea_m2, order.by = DTM)

  don <- cbind(ref,measure)
  colnames(don) <- c("Odtok bez opatreni", "Odtok s opatrenim")
  mmain=paste0("Graf simulace odtoku z povodi ")
  p <- dygraph(don,main = mmain) %>%
    dyAxis("y", label = "Prumerny denni odtok v m3/den") %>%
    dyOptions(labelsUTC = TRUE, fillGraph=TRUE, fillAlpha=0.1, drawGrid = FALSE, colors = c("BLACK","RED","BLUE","GREEN")) %>%
    dyRangeSelector() %>%
    dyCrosshair(direction = "vertical") %>%
    dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE)  %>%
    dyRoller(rollPeriod = 1)
  if (uloz_graf){
    saveWidget(p, file=paste0( getwd(), "/",mmain,".html"))
  }
  p
}



graf_transpirace_WTLND<-function(outDF_sim_measure,CatchmentArea_m2,uloz_graf){

  DTM=as.Date(paste(outDF_sim_measure$YEAR, outDF_sim_measure$MONTH, outDF_sim_measure$DAY, sep = "-"))
  Etpond  = outDF_sim_measure$TRNS/1000*CatchmentArea_m2

  Vypar <- xts(x = Etpond, order.by = DTM)
  don <- cbind(Vypar)
  colnames(don) <- c("Vypar z opatreni")

  mmain=paste0("Graf vyparu z opatreni ")
  p <- dygraph(don,main = mmain) %>%
    dyAxis("y", label = "Denní vypar v m3/den") %>%
    dyOptions(labelsUTC = TRUE, fillGraph=TRUE, fillAlpha=0.1, drawGrid = FALSE, colors = c("BLACK","RED","BLUE","GREEN")) %>%
    dyRangeSelector() %>%
    dyCrosshair(direction = "vertical") %>%
    dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE)  %>%
    dyRoller(rollPeriod = 1)
  if (uloz_graf){
    saveWidget(p, file=paste0( getwd(), "/",mmain,".html"))
  }
  p
}

graf_porovnani_celkoveho_vyparu<-function(outDF_ref,outDF_sim_measure,uloz_graf){

  DTM=as.Date(paste(outDF_sim_measure$YEAR, outDF_sim_measure$MONTH, outDF_sim_measure$DAY, sep = "-"))
  AET_measure  = outDF_sim_measure$AET
  AET_NOmeasure  = outDF_ref$AET

  Vypar_s <- xts(x = AET_measure, order.by = DTM)
  Vypar_bez <- xts(x = AET_NOmeasure, order.by = DTM)

  don <- cbind(Vypar_s,Vypar_bez)
  colnames(don) <- c("Vypar s opatrenim", "Vypar bez opatreni")

  mmain=paste0("Porovnani vyparu z celeho povodi ")
  p <- dygraph(don,main = mmain) %>%
    dyAxis("y", label = "Denní výpar v mm/den") %>%
    dyOptions(labelsUTC = TRUE, fillGraph=TRUE, fillAlpha=0.1, drawGrid = FALSE, colors = c("BLACK","RED","BLUE","GREEN")) %>%
    dyRangeSelector() %>%
    dyCrosshair(direction = "vertical") %>%
    dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE)  %>%
    dyRoller(rollPeriod = 1)
  if (uloz_graf){
    saveWidget(p, file=paste0( getwd(), "/",mmain,".html"))
  }
  p
}


graf_porovnani_vyparu_z_retence<-function(outDF_ref,outDF_sim_measure,uloz_graf){

  DTM=as.Date(paste(outDF_sim_measure$YEAR, outDF_sim_measure$MONTH, outDF_sim_measure$DAY, sep = "-"))
  measure  = outDF_sim_measure$ETSW+outDF_sim_measure$TRNS
  NOmeasure  = outDF_ref$ETSW

  Vypar_s <- xts(x = measure, order.by = DTM)
  Vypar_bez <- xts(x = NOmeasure, order.by = DTM)

  don <- cbind(Vypar_s,Vypar_bez)
  colnames(don) <- c("Vypar s opatrenim", "Vypar bez opatreni")

  mmain=paste0("Porovnani vyparu z celeho povodi ")
  p <- dygraph(don,main = mmain) %>%
    dyAxis("y", label = "Denní výpar v mm/den") %>%
    dyOptions(labelsUTC = TRUE, fillGraph=TRUE, fillAlpha=0.1, drawGrid = FALSE, colors = c("BLACK","RED","BLUE","GREEN")) %>%
    dyRangeSelector() %>%
    dyCrosshair(direction = "vertical") %>%
    dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE)  %>%
    dyRoller(rollPeriod = 1)
  if (uloz_graf){
    saveWidget(p, file=paste0( getwd(), "/",mmain,".html"))
  }
  p
}
