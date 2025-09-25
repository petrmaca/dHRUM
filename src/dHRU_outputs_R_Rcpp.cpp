#include <Rcpp.h>
#include "dHRUM.h"

//' Provides dHRUM outputs - time series
//'
//' return matrix with state variables and fluxes averaged over basin area.
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @return list with matrix of caldata \code{[,1:4]} and hdata ts variables \code{[,5:27]} and names of vars
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
//' calcHBInAlldHrus(dHRUM_ptr = dhrus)
//' gatherHBdata(dHRUM_ptr = dhrus)
//' outDta <- getOutput(dHRUM_ptr = dhrus)
// [[Rcpp::export]]
Rcpp::List getOutput(Rcpp::XPtr<dHRUM> dHRUM_ptr){
  unsigned nrowOutMat = dHRUM_ptr.get()->get_numTS();
  unsigned ncolOutMat = numTSvars+4;
  Rcpp::NumericMatrix outDta( nrowOutMat, ncolOutMat ) ;
  hdata helpVal = dHRUM_ptr.get()->get_HbDta(all_ts[0]);
  // numTSvar
#pragma omp parallel for
  for(unsigned j=0;j<4;j++){
    caldata helpVal = dHRUM_ptr.get()->get_CalDta(all_caDT[j]);
    for(unsigned i=0; i<nrowOutMat; i++){
      outDta(i,j) = helpVal[i] ;
    }
  }

#pragma omp parallel for
  for(unsigned j=4;j<ncolOutMat;j++){
    hdata helpVal = dHRUM_ptr.get()->get_HbDta(all_ts[j-4]);
    for(unsigned i=0; i<nrowOutMat; i++){
      outDta(i,j) = helpVal[i] ;
    }
  }
  Rcpp::StringVector VarsNams({"YEAR", "MONTH", "DAY", "JDAY",    \
                              "PREC","SNOW","AET","PET","TEMP",   \
                              "MELT","TROF","STEF","CANF","CANS", \
                              "STES","EVAC","EVAS","EVBS","INTS", \
                              "SOIS","GROS","SURS","TOTR","BASF",\
                              "DIRR","PERC","PREF", "ETSW","PONS"});

  return Rcpp::List::create(
    Rcpp::Named("outDta") = outDta,
    Rcpp::Named("VarsNams") = VarsNams
  );
}

//' Write Output of overall dHru simulation to file
//'
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @param namOutFilet file and path to write output to
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dhrus,as.numeric(ParDF[1,]),names(ParDF))
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
//' calcHBInAlldHrus(dHRUM_ptr = dhrus)
//' gatherHBdata(dHRUM_ptr = dhrus)
//' #printToFile(dHRUM_ptr = dhrus, "path/filename")
// [[Rcpp::export]]
void printToFile(Rcpp::XPtr<dHRUM> dHRUM_ptr, std::string namOutFilet) {

  return dHRUM_ptr.get()->printAllDta(namOutFilet);

}

//' Provides dHRUM outputs - time series for each HRU
//'
//' return list matrix with state variables and fluxes averaged over basin area.
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @return list with matrix of caldata \code{[,1:4]} and hdata ts variables \code{[,5:27]} and names of vars
//' @export
//' @examples
//' nHrus <- 2
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'                       SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
//' calcHBInAlldHrus(dHRUM_ptr = dhrus)
//' gatherHBdata(dHRUM_ptr = dhrus)
//' outDta <- getOutputDist(dHRUM_ptr = dhrus)
//' outDta <- getOutputDist(dHRUM_ptr = dhrus)
//' outDF <- cbind(outDta$outDta, outDta$Ids)
//' outDF <- data.frame(outDF)
//' names(outDF) <-c(outDta$VarsNams,"HruIDs")
// [[Rcpp::export]]
Rcpp::List getOutputDist(Rcpp::XPtr<dHRUM> dHRUM_ptr){
  unsigned nrowOutMat = dHRUM_ptr.get()->get_numTS() * dHRUM_ptr->getdHRUdim();
  unsigned ncolOutMat = numTSvars+4;
  Rcpp::NumericMatrix outDta( nrowOutMat, ncolOutMat ) ;
  Rcpp::StringVector Ids(nrowOutMat);

  hdata helpVal = dHRUM_ptr.get()->get_HbDta(all_ts[0]);

  // building the vector with IDs
  unsigned indexrow=0;
  for(unsigned itHRU=0; itHRU<dHRUM_ptr->getdHRUdim(); itHRU++){
    for(unsigned i=0;i<dHRUM_ptr.get()->get_numTS();i++){
      Ids[indexrow] = dHRUM_ptr.get()->getSingleHruId(itHRU);
      indexrow++;
      }
    }
  // numTSvar

#pragma omp parallel for
  for(unsigned j=0;j<4;j++){
    unsigned indexrowCal=0;
    for(unsigned itHRU=0; itHRU<dHRUM_ptr->getdHRUdim(); itHRU++){
      caldata helpVal = dHRUM_ptr.get()->get_CalDta(all_caDT[j]);
      for(unsigned i=0; i<dHRUM_ptr.get()->get_numTS(); i++){
        outDta(indexrowCal,j) = helpVal[i];
        indexrowCal++;
      }
    }
  }

#pragma omp parallel for
  for(unsigned j=4;j<ncolOutMat;j++){
    unsigned indexrowNu=0;
    for(unsigned itHRU=0; itHRU<dHRUM_ptr->getdHRUdim(); itHRU++){
      hdata helpVal = dHRUM_ptr.get()->get_HbDta(all_ts[j-4]);
      for(unsigned i=0; i<dHRUM_ptr.get()->get_numTS(); i++){
        outDta(indexrowNu,j) = helpVal[i] ;
        indexrowNu++;
      }
    }
  }

  Rcpp::StringVector VarsNams({"YEAR", "MONTH", "DAY", "JDAY",    \
                              "PREC","SNOW","AET","PET","TEMP",   \
                              "MELT","TROF","STEF","CANF","CANS", \
                              "STES","EVAC","EVAS","EVBS","INTS", \
                              "SOIS","GROS","SURS","TOTR","BASF",\
                              "DIRR","PERC","PREF", "ETSW","PONS"});

  return Rcpp::List::create(
    Rcpp::Named("outDta") = outDta,
    Rcpp::Named("VarsNams") = VarsNams,
    Rcpp::Named("Ids") = Ids
  );
}


//' Run dHRUM and provides with model outputs
//'
//' return list matrix with state variables and fluxes averaged over basin area.
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @return list with matrix of caldata \code{[,1:4]} and hdata ts variables \code{[,5:27]} and names of vars
//' @export
//' @examples
//' nHrus <- 2
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'                       SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
//' outDta <- dHRUMrun(dHRUM_ptr = dhrus)
//' outDF <- data.frame(outDta$outDta)
//' names(outDF) <-c(outDta$VarsNams)
// [[Rcpp::export]]
Rcpp::List dHRUMrun(Rcpp::XPtr<dHRUM> dHRUM_ptr){

  dHRUM_ptr.get()->calcHbToAllHrus();
  dHRUM_ptr.get()->gatherTsFromHrus();

  unsigned nrowOutMat = dHRUM_ptr.get()->get_numTS();
  unsigned ncolOutMat = numTSvars+4;
  Rcpp::NumericMatrix outDta( nrowOutMat, ncolOutMat ) ;
  hdata helpVal = dHRUM_ptr.get()->get_HbDta(all_ts[0]);
  // numTSvar
  for(unsigned j=0;j<4;j++){
    caldata helpVal = dHRUM_ptr.get()->get_CalDta(all_caDT[j]);
    for(unsigned i=0; i<nrowOutMat; i++){
      outDta(i,j) = helpVal[i] ;
    }
  }

  for(unsigned j=4;j<ncolOutMat;j++){
    hdata helpVal = dHRUM_ptr.get()->get_HbDta(all_ts[j-4]);
    for(unsigned i=0; i<nrowOutMat; i++){
      outDta(i,j) = helpVal[i] ;
    }
  }
  Rcpp::StringVector VarsNams({"YEAR", "MONTH", "DAY", "JDAY",    \
                              "PREC","SNOW","AET","PET","TEMP",   \
                              "MELT","TROF","STEF","CANF","CANS", \
                              "STES","EVAC","EVAS","EVBS","INTS", \
                              "SOIS","GROS","SURS","TOTR","BASF", \
                              "DIRR","PERC","PREF", "ETWS","PONS"});
  return Rcpp::List::create(
    Rcpp::Named("outDta") = outDta,
    Rcpp::Named("VarsNams") = VarsNams
  );
}


//' Run dHRUM and provides with distributed model outputs obtained form each HRU
//'
//' return list matrix with state variables and fluxes from each HRU.
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @return list with matrix of caldata \code{[,1:4]} and hdata ts variables \code{[,5:27]} and names of vars
//' @export
//' @examples
//' nHrus <- 2
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTInputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'                       SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
//' outDta <- dHRUMrunDist(dHRUM_ptr = dhrus)
//' outDF <- data.frame(outDta$outDta)
//' outDF <- cbind(outDF, outDta$Ids)
//' names(outDF) <-c(outDta$VarsNams)
// [[Rcpp::export]]
Rcpp::List dHRUMrunDist(Rcpp::XPtr<dHRUM> dHRUM_ptr){

  dHRUM_ptr.get()->calcHbToAllHrus();
  dHRUM_ptr.get()->gatherTsFromHrus();
  unsigned nrowOutMat = dHRUM_ptr.get()->get_numTS() * dHRUM_ptr->getdHRUdim();
  unsigned ncolOutMat = numTSvars+4;
  Rcpp::NumericMatrix outDta( nrowOutMat, ncolOutMat ) ;
  Rcpp::StringVector Ids(nrowOutMat);

  hdata helpVal = dHRUM_ptr.get()->get_HbDta(all_ts[0]);

  unsigned indexrow=0;
  for(unsigned itHRU=0; itHRU<dHRUM_ptr->getdHRUdim(); itHRU++){
    for(unsigned i=0;i<dHRUM_ptr.get()->get_numTS();i++){
      Ids[indexrow] = dHRUM_ptr.get()->getSingleHruId(itHRU);
      indexrow++;
    }
  }
  // numTSvar

  for(unsigned j=0;j<4;j++){
    unsigned indexrowCal=0;
    for(unsigned itHRU=0; itHRU<dHRUM_ptr->getdHRUdim(); itHRU++){
      caldata helpVal = dHRUM_ptr.get()->get_CalDta(all_caDT[j]);
      for(unsigned i=0; i<dHRUM_ptr.get()->get_numTS(); i++){
        outDta(indexrowCal,j) = helpVal[i];
        indexrowCal++;
      }
    }
  }

  for(unsigned j=4;j<ncolOutMat;j++){
    unsigned indexrowNu=0;
    for(unsigned itHRU=0; itHRU<dHRUM_ptr->getdHRUdim(); itHRU++){
      for(unsigned i=0; i<dHRUM_ptr.get()->get_numTS(); i++){
        numberSel helpVal = dHRUM_ptr.get()->getTsDta(all_ts[j-4], itHRU, i);
        outDta(indexrowNu,j) = helpVal;
        indexrowNu++;
      }
    }
  }

  Rcpp::StringVector VarsNams({"YEAR", "MONTH", "DAY", "JDAY",    \
                              "PREC","SNOW","AET","PET","TEMP",   \
                              "MELT","TROF","STEF","CANF","CANS", \
                              "STES","EVAC","EVAS","EVBS","INTS", \
                              "SOIS","GROS","SURS","TOTR","BASF",\
                              "DIRR","PERC","PREF", "ETWS","PONS","HruIds"});

  return Rcpp::List::create(
    Rcpp::Named("outDta") = outDta,
    Rcpp::Named("VarsNams") = VarsNams,
    Rcpp::Named("Ids") = Ids
  );
}
