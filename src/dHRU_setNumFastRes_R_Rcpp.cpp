#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the similar values of params to dHRU model for all single HRUs.
//'
//' Setting the groudnwater type of dHRUM equal to all HRUs. Possibe types: \code{LIN_RES} \code{LINL_RES} \code{LINBY_RES}
//' \code{POW_RES}, \code{EXP_RES} \code{LIN_2SE} \code{LIN_2PA} \code{FLEX_RES}
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param gwTypes a charater vector of GW type names
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' ups = 1:nHrus
//' setNumFastResAlldHrus(dHRUM_ptr = dhrus,numFastRes=ups,hruIds=IdsHrus)
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
