#include <Rcpp.h>
#include "dHRUM.h"

//' Calculates the values of fluxes and state variables for all Hrus in dHRUM
//'
//' Updating the states ad values of all ts variables in all single Hrus of dHRUM
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' filname2 = "../../PDM/Development/PDM_dist/data/tests/inALL/BP_1960_01_01.txt"
//' setInputsToAlldHrus(dHRUM_ptr = dhrus, filname2)
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF),TRUE,0)
//' calcPetToHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
//' calcHBInAlldHrus(dHRUM_ptr = dhrus)
// [[Rcpp::export]]
void calcHBInAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr) {

  return dHRUM_ptr.get()->calcHbToAllHrus();

}

//' Calculates catchment spatially averaged hydrological balance
//'
//' Calculates catchment spatial average for fluxes and state variables
//' after the HB on every single HRU is updated.
//'
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' filname2 = "../../PDM/Development/PDM_dist/data/tests/inALL/BP_1960_01_01.txt"
//' setInputsToAlldHrus(dHRUM_ptr = dhrus, filname2)
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF),TRUE,0)
//' calcPetToHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
//' calcHBInAlldHrus(dHRUM_ptr = dhrus)
//' gatherHBdata(dHRUM_ptr = dhrus)
// [[Rcpp::export]]
void gatherHBdata(Rcpp::XPtr<dHRUM> dHRUM_ptr) {

  return dHRUM_ptr.get()->gatherTsFromHrus();

}
