
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


graf_porovnani_TOTR<-function(outDF_sim_pond,PondArea_m2,CatchmentArea_m2,uloz_graf){
  DTM=as.Date(paste(outDF_sim_pond$YEAR, outDF_sim_pond$MONTH, outDF_sim_pond$DAY, sep = "-"))
  
  no_pond <- xts(x = (outDF_sim_pond$BASF+outDF_sim_pond$DIRR)/1000*CatchmentArea_m2, order.by = DTM)
  pond <- xts(x = outDF_sim_pond$TOTR/1000*CatchmentArea_m2, order.by = DTM)
  
  don <- cbind(no_pond,pond)
  colnames(don) <- c("Odtok bez nádrže", "Odtok z nádrže")
  mmain=paste0("Graf simulace odtoku z povodi ")
  p <- dygraph(don,main = mmain) %>%
    dyAxis("y", label = "Průměrný denní odtok v m3/den") %>%
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


graf_porovnani_pritoku_zasoby<-function(outDF_sim_pond,PondArea_m2,CatchmentArea_m2,uloz_graf){
  pritok_m3_den=((outDF_sim_pond$BASF+outDF_sim_pond$DIRR)/1000*CatchmentArea_m2)*-1 #pritok invertuju, aby mohl jit zhoza dolu v grfu
  objem_3m=outDF_sim_pond$PONS
  DTM=as.Date(paste(outDF_sim_pond$YEAR, outDF_sim_pond$MONTH, outDF_sim_pond$DAY, sep = "-"))
  
  Pritok <- xts(x = pritok_m3_den, order.by = DTM)
  Zasoba <- xts(x = objem_3m, order.by = DTM)

  
  don <- cbind(Zasoba,Pritok )
  colnames(don) <- c("Objem nádrže", "Přítok")
  mmain=paste0("Bilanční graf nádrže (Objem vs. Přítok)")
  p <-dygraph(don, main = mmain) %>%
    # Objem nádrže bude na levé ose (Y)
    dySeries("Objem nádrže", axis = 'y', color = "darkblue", strokeWidth = 2) %>%
    # Přítok bude na pravé ose (Y2) a vykreslí se jako výplň (bary/plocha)
    dySeries("Přítok", axis = 'y2', color = "rgba(0, 150, 255, 0.6)", stepPlot = TRUE, fillGraph = TRUE) %>%
    # Nastavení levé osy (Objem) - klasická osa od 0 nahoru
    dyAxis("y", label = "Objem vody v nádrži [m³]") %>%
    # Nastavení pravé osy (Přítok) - klíčový trik s 'valueFormatter' a 'axisLabelFormatter'
    dyAxis("y2", 
           label = "Přítok [m³/s]", 
           independentTicks = TRUE,
           # Tento JavaScript kód schová záporné znaménko na ose
           axisLabelFormatter = "function(v) { return Math.abs(v); }",
           # Tento JavaScript kód schová záporné znaménko v legendě při najetí myší
           valueFormatter = "function(v) { return Math.abs(v).toFixed(2); }") %>%
    # Přidání spodního posuvníku pro zoomování
    dyRangeSelector(height = 40) %>%
    # Zlepšení vzhledu legendy
    dyLegend(width = 400, show = "always", hideOnMouseOut = FALSE)
  
  if (uloz_graf){
    saveWidget(p, file=paste0( getwd(), "/",mmain,".html"))  
  }
  p
}


graf_vypar_z_nadrze<-function(outDF_sim_pond,PondArea_m2,uloz_graf){
  
  DTM=as.Date(paste(outDF_sim_pond$YEAR, outDF_sim_pond$MONTH, outDF_sim_pond$DAY, sep = "-"))
  Etpond  = outDF_sim_pond$ETPO/1000*PondArea_m2
  
  Vypar <- xts(x = Etpond, order.by = DTM)
  don <- cbind(Vypar)
  colnames(don) <- c("Výpar z nádrže")
  
  mmain=paste0("Graf výparu z nádrže ")
  p <- dygraph(don,main = mmain) %>%
    dyAxis("y", label = "Denní výpar v m3/den") %>%
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

graf_porovnani_vyparu<-function(outDF_sim,outDF_sim_pond,PondArea_m2,CatchmentArea_m2=Area,uloz_graf){
  
  DTM=as.Date(paste(outDF_sim_pond$YEAR, outDF_sim_pond$MONTH, outDF_sim_pond$DAY, sep = "-"))
  AET_pond  = outDF_sim_pond$AET
  AET_NOpond  = outDF_sim$AET
  
  Vypar_s <- xts(x = AET_pond, order.by = DTM)
  Vypar_bez <- xts(x = AET_NOpond, order.by = DTM)
  
  don <- cbind(Vypar_s,Vypar_bez)
  colnames(don) <- c("Výpar s nádrží", "Výpar bez nádrže")
  
  mmain=paste0("Porovnání výparů z celého povodí ")
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