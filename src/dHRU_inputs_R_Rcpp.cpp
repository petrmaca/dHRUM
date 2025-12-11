#include <Rcpp.h>
#include "dHRUM.h"

//' Sets similar input data obtained from file to all single Hrus at dHRUM instance.
//'
//' Loads the data from file to a single dHRUM instance created  the \code{initdHruModel(nHrus,Areas,IdsHrus)} function
//' All inputs are same for each single HRU unit. File has on its first row YYYY MM DD,
//' remaining columns should have Temperature and Precipitation data.
//'
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @param namInpFilet a character vector to a single file with input data to dHRUM.
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus,5)
//' filname2 = "../Calibrations/Amalie/indata/BP_1960_01_01.txt"
//' setPTInputsToAlldHrusFromFile(dHRUM_ptr = dhrus, filname2)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
void setPTInputsToAlldHrusFromFile(Rcpp::XPtr<dHRUM> dHRUM_ptr, std::string namInpFilet) {
  // Rcpp::Rcout << "Input data fully loaded.";
  // See comment on returning pointer
  //https://stackoverflow.com/questions/59384221/proper-way-to-return-a-pointer-to-a-new-object-from-an-rcpp-function
  //ToDo rcpp modules and factory wrap of pointer
  return dHRUM_ptr.get()->setInputsToAllHrus(namInpFilet);

}

//' Sets the similar Precipitation, Temperature vectors to dHRUM and init's the date using beg. of period.
//'
//' Setting the similar vector of Precipitation and temperature to all single HRU.
//' Setting the calender using the first date fo period using the first date of period
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param Prec vector of values of precipitation
//' @param Temp vector of values of temperature
//' @param inDate the first date of simulation period
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
// [[Rcpp::export]]
void setPTInputsToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector Prec, Rcpp::NumericVector Temp, Rcpp::Date inDate) {

  unsigned Myear = 0, Mmonth = 0, Mday = 0;

  Myear = (unsigned) inDate.getYear();
  Mmonth = (unsigned) inDate.getMonth();
  Mday = (unsigned) inDate.getDay();
  // Rcpp::Rcout << Myear << " mm " << Mmonth << " dd " << Mday+1 <<"\n";
  unsigned ndatPrec = 0, ndatTemp = 0;
  ndatPrec = Prec.size();
  ndatTemp = Temp.size();

  if((ndatTemp!=ndatPrec)) {
    Rcpp::Rcout << "The size of Temp vector is " << ndatTemp <<"  and size of Precip is "<< ndatPrec << "\n";
      Rcpp::stop("\n Different size of input precipitation and temperature data.\n");
  }

  hdata mPrec(1,1), mTemp(1,1);

  mPrec.resize(ndatPrec);
  mTemp.resize(ndatPrec);

  for(unsigned it=0;it<ndatPrec;it++){
    mPrec[it] = Prec[it];
    mTemp[it] = Temp[it];
  }
  // for(unsigned it=0;it<ndatPrec;it++){
  //   Rcpp::Rcout << mPrec[it] << "\n";
  //   Rcpp::Rcout << mTemp[it] << "\n";
  // }
  // Rcpp::Rcout <<  "Loading the gw type  to HRU ID "  << std::endl;
  //    dHruVec[it].set_paramsToSim(parsToLoad);
  // Rcpp::Rcout<<"threads="<<omp_get_num_threads()<<std::endl;

  dHRUM_ptr.get()->loadPTDatToAllHrus(mPrec, mTemp, 0.0, Myear, Mmonth, Mday);

  return  ;
}

//' Sets the similar Precipitation, Temperature, and Date vectors to dHRU.
//'
//' Setting the similar vectors of Precipitation and Temperature
//' to all single HRUs of dHRU. Loading the Date vector  from vector like variable.
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @param Prec vector of values of precipitation
//' @param Temp vector of values of temperature
//' @param DateVec vector of dates \code{as.Date("YYYY-MM-DD")}
//' @export
//' @examples
//' nHrus <- 2
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus,2)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
// [[Rcpp::export]]
void setPTDateInputsToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector Prec, Rcpp::NumericVector Temp, Rcpp::DateVector DateVec) {

  // unsigned Myear = 0, Mmonth = 0, Mday = 0;
  //
  // Myear = (unsigned) inDate.getYear();
  // Mmonth = (unsigned) inDate.getMonth();
  // Mday = (unsigned) inDate.getDay();
  // Rcpp::Rcout << Myear << " mm " << Mmonth << " dd " << Mday+1 <<"\n";
  unsigned ndatPrec = 0, ndatTemp = 0, ndatDate = 0;
  ndatPrec = Prec.size();
  ndatTemp = Temp.size();
  ndatDate = DateVec.size();

  if((ndatTemp!=ndatPrec)&&(ndatPrec!=ndatDate)) {
    Rcpp::stop("\n Different size of input precipitation, temperature, and dates vectors.\n");
  }

  hdata mPrec(1,1), mTemp(1,1);
  caldata myear(1,1), mmonth(1,1), mday(1,1);

  mPrec.resize(ndatPrec);
  mTemp.resize(ndatPrec);
  myear.resize(ndatPrec);
  mmonth.resize(ndatPrec);
  mday.resize(ndatPrec);

  for(unsigned it=0;it<ndatPrec;it++){
    mPrec[it] = Prec[it];
    mTemp[it] = Temp[it];
    // myear[it] = (unsigned) DateVec[it].getYear();
    Rcpp::Date myDat = Rcpp::Date(DateVec[it]);
    myear[it] = (unsigned) myDat.getYear();
    mmonth[it] = (unsigned) myDat.getMonth();
    mday[it] = (unsigned) myDat.getDay();
  }

  //   for(unsigned it=0;it<ndatPrec;it++){
  //     Rcpp::Rcout << mmonth[it] << "\n";
  //     Rcpp::Rcout << mmonth[it] << "\n";
  //   }
  dHRUM_ptr.get()->load_PrecTempToAllHrus(mPrec, mTemp);
  dHRUM_ptr.get()->load_CalDataToAllHrus(myear, mmonth, mday);
  dHRUM_ptr.get()->initdHRUbasinDTA();

  return  ;

}

//' Sets the distributed Precipitation, Temperature vectors to distributed dHRUM and init's the date using beg. of period.
//'
//' Setting the different vector of Precipitation and temperature to all single HRU.
//' Setting the calender using the first date fo period using the first date of period
//' The ordering og Ids must be constant for all input data uploads
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param DataDF dataframe with DTM, Precipitation, Temperature, and HRU Ids
//' @export
//' @examples
//' nHrus <- 2
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus,2)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
void setPTInputsToDistdHRUM(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::DataFrame DataDF) {

  Rcpp::CharacterVector HruIdVec = DataDF["HruId"];
  Rcpp::NumericVector Prec = DataDF["P"];
  Rcpp::NumericVector Temp = DataDF["T"];
  Rcpp::DateVector DateVec = DataDF["DTM"];

  unsigned nDatInOneHru = 0, nHrusDF = 0, ndat =0;

  std::unordered_set<SEXP> uniqueHRUs(HruIdVec.begin(),HruIdVec.end());
  nHrusDF = uniqueHRUs.size();

  // Rcpp::Rcout << "the number of HruId in Df: " << nHrusDF << "\n";
  if(nHrusDF != ((unsigned) dHRUM_ptr.get()->getdHRUdim())){
    Rcpp::stop("\n Different number of Hru's in data.framne and dHRUM.\n");
  }

  nDatInOneHru = std::count(HruIdVec.begin(),HruIdVec.end(),HruIdVec[0]);
  // std::cout << ndat;
  ndat = nDatInOneHru * nHrusDF;

  if(ndat != ((unsigned) DataDF.nrows())){
    Rcpp::stop("\n Different and non constant number of ts data in data.framne for Hrus.\n");
  }
  // Rcpp::Rcout << "The number of all ts data: " << ndat << "\n";

  Rcpp::Date myDat = Rcpp::Date(DateVec[0]);

  unsigned myear= 0, mmonth = 0, mday = 0;
  myear = (unsigned) myDat.getYear();
  mmonth = (unsigned) myDat.getMonth();
  mday = (unsigned) myDat.getDay();

  unsigned TScounter = 0;
  for(unsigned int itHru=0; itHru < nHrusDF;itHru++) {
    hdata mPrec(1,1), mTemp(1,1);
    mPrec.resize(nDatInOneHru);
    mTemp.resize(nDatInOneHru);
    for(unsigned int ts=0;ts<nDatInOneHru;ts++) {
      mPrec[ts] = Prec[TScounter];
      mTemp[ts] = Temp[TScounter];
      // myear[it] = (unsigned) DateVec[it].getYear();
      TScounter++;
    }
    dHRUM_ptr.get()->loadPTInputsToOneHru(mPrec,mTemp,0,myear,mmonth,mday,itHru);
  }

  dHRUM_ptr.get()->initdHRUbasinDTA();

  return ;

}
