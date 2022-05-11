#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the number of fast runoff reservoirs to dHRU model for all single HRUs.
//'
//' Setting the number of reservoir for direct runoff on dHRUM  HRUs.
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param numFastRes vector with integer numbers
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' ups = 1:nHrus
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setNumFastResAlldHrus(dHRUM_ptr = dhrus,numFastRes=ups,hruIds=IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
void setNumFastResAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector numFastRes, Rcpp::CharacterVector hruIds) {
  unsigned dimdHRUM = 0;

  unsigned numNFR = numFastRes.size();
  unsigned numHruIdNames = hruIds.size();

  dimdHRUM = dHRUM_ptr.get()->getdHRUdim();


  if(dimdHRUM!=numHruIdNames) {
    Rcpp::Rcout << "The dim of dHRUs is different than number of Ids provided   " << numHruIdNames <<"\n";
    Rcpp::stop("\nWrong size of numbers Ids and dim DHRUM.\n");
  } else if(numNFR!=numHruIdNames) {
    Rcpp::Rcout << "The dim of vector with numbers of fast reservoirs does not correspond to the number of HRUs " << numNFR <<"\n";
    Rcpp::stop("\nWrong size of numbers of fast reservoirs types.\n");
  } else {
    caldata helpVec(1,1);
    helpVec.resize(numNFR);
    for(unsigned it=0;it<numNFR;it++){
      helpVec[it] = numFastRes[it];
    }
    dHRUM_ptr.get()->set_numFastReservoirs(helpVec);
    dHRUM_ptr.get()->set_numFastReservoirsToHrus();
  }

  return;
}
