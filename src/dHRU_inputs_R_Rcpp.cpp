#include <Rcpp.h>
#include "dHRUM.h"

//' Sets similar input data obtained from file to all single Hrus at dHRUM instance.
//'
//' Loads the data from file to a single dHRU instance created  the \code{initdHruModel(nHrus,Areas,IdsHrus)} function
//' All iputs are same for each single HRU unit. File has on its first row YYYY MM DD,
//' remaining columns should have Temperature and Precipitation data.
//'
//'
//' @param dHRU_ptr pointer to dHRU instance
//' @param namInpFilet a chacter vector to a single file with input data to dHRUM.
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' filname2 = "../../PDM/Development/PDM_dist/data/tests/inALL/BP_1960_01_01.txt"
//' setInputsToAlldHrus(dHRUM_ptr = dhrus, filname2)
// [[Rcpp::export]]
void setInputsToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, std::string namInpFilet) {
  // Rcpp::Rcout << "Input data fully loaded.";
  // See comment on returning pointer
  //https://stackoverflow.com/questions/59384221/proper-way-to-return-a-pointer-to-a-new-object-from-an-rcpp-function
  //ToDo rcpp modules and factory wrap of pointer
  return dHRUM_ptr.get()->setInputsToAllHrus(namInpFilet);

}

//' Sets the Precipitation, Temperature vectors to dHRUM and init's the date using beg. of period.
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
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTinputsToAlldHrus(dhrus, Prec = prec, Temp = temp, as.Date("1990/01/30"))
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
  dHRUM_ptr.get()->loadPTDatToAllHrus(mPrec, mTemp, 0.0, Myear, Mmonth, Mday);

  return  ;
}

//' Sets the Precipitation, Temperature, and Date vectors to dHRU.
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
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
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