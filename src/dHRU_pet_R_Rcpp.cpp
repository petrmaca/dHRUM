#include <Rcpp.h>
#include "dHRUM.h"

//' Calculates the values of Potential evapotranspiration on all singleHrus
//'
//' Setting of Pet method is done by \code{PetTypeStr}, methods implemented are: \code{Oudin}, \code{Hamon}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param Latitude single number for Oudin method
//' @param PetTypeStr variable on selection of PET models
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' filname2 = "../../PDM/Development/PDM_dist/data/tests/inALL/BP_1960_01_01.txt"
//' setInputsToAlldHrus(dHRUM_ptr = dhrus, filname2,TRUE,0)
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF),TRUE,0)
//' calcPetToHrus(dHRUM_ptr = dhrus,50.1,"Hamon",TRUE,0)
// [[Rcpp::export]]
void calcPetToAllHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, numberSel Latitude, std::string PetTypeStr){
  std::vector<std::string> PettypesVec {"Oudin", "Hamon"};
  bool petBB=false;
  for(unsigned pp=0;pp<2;pp++){
    if(PettypesVec[pp] == PetTypeStr){
      petBB = true;
    }
  }
  if(!petBB) {
    Rcpp::stop("\n Wrong calling to the Pet Type method, try 'Oudin' or 'Hamon'.\n");
  }

  std::map<std::string, pet_Type> s_mapStringToPet_Type = {
    {"Oudin", pet_Type::OUDIN},
    {"Hamon", pet_Type::HAMON}
  };
  pet_Type  myPetType;
  switch(s_mapStringToPet_Type[PetTypeStr]) {
  case pet_Type::OUDIN:
    myPetType = pet_Type::OUDIN;
    // Rcpp::Rcout << "\n oudinddd \n";
    break;
  case pet_Type::HAMON:
    myPetType = pet_Type::HAMON;
    // Rcpp::Rcout << "\n hamon \n";
    break;
  }

  return dHRUM_ptr.get()->calcPetToAllHrus(Latitude, myPetType);

}
