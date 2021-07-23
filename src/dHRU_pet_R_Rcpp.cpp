#include <Rcpp.h>
#include "dHRUM.h"

//' Calculates the values of Potential evapotranspiration on all singleHrus taking constant latitude
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
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
// [[Rcpp::export]]
void calcPetToAllHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, numberSel Latitude, std::string PetTypeStr){
  std::vector<std::string> PettypesVec {"Oudin", "Hamon","Thornthwaite","BlaneyCriddle","Jensenhaise","McQuinnessbordne"};
  bool petBB=false;
  for(unsigned pp=0;pp<6;pp++){
    if(PettypesVec[pp] == PetTypeStr){
      petBB = true;
    }
  }
  if(!petBB) {
    Rcpp::stop("\n Wrong calling to the Pet Type method, try 'Oudin' or 'Hamon or 'Thornthwaite' or 'BlaneyCriddle' or 'Jensenhaise' or 'McQuinnessbordne'.\n");
  }

  std::map<std::string, pet_Type> s_mapStringToPet_Type = {
    {"Oudin", pet_Type::OUDIN},
    {"Hamon", pet_Type::HAMON},
    {"Thornthwaite", pet_Type::THORNTHWAITE},
    {"BlaneyCriddle", pet_Type::BLANEYCRIDDLE},
    {"Jensenhaise", pet_Type::JENSENHAISE},
    {"McQuinnessbordne", pet_Type::MCGUINNESSBORDNE}
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
  case pet_Type::THORNTHWAITE:
    myPetType = pet_Type::THORNTHWAITE;
    // Rcpp::Rcout << "\n hamon \n";
    break;
  case pet_Type::BLANEYCRIDDLE:
    myPetType = pet_Type::BLANEYCRIDDLE;
    // Rcpp::Rcout << "\n hamon \n";
    break;
  case pet_Type::JENSENHAISE:
    myPetType = pet_Type::JENSENHAISE;
    // Rcpp::Rcout << "\n hamon \n";
    break;
  case pet_Type::MCGUINNESSBORDNE:
    myPetType = pet_Type::MCGUINNESSBORDNE;
    // Rcpp::Rcout << "\n hamon \n";
    break;
  }

  return dHRUM_ptr.get()->calcPetToAllHrus(Latitude, myPetType);

}

//' Calculates the values of Potential evapotranspiration on all singleHrus taking nonconstant latitude
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
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
// [[Rcpp::export]]
void calcPetToAllHrusDist(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector Latitude, Rcpp::CharacterVector PetTypeStrNames, Rcpp::CharacterVector HruIds){
  std::vector<std::string>  petNameStr = Rcpp::as<std::vector<std::string> >(PetTypeStrNames);
  unsigned dHRUMdim = dHRUM_ptr->getdHRUdim();

  for(unsigned it=0; it<dHRUMdim;it++ ){
    if ( std::find(allPetNames.begin(), allPetNames.end(), petNameStr[it]) == allPetNames.end()) {
      Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
      Rcpp::stop("\n Wrong names of Pet names.\n");
    }
  }


//
//
//   std::map<std::string, pet_Type> s_mapStringToPet_Type = {
//     {"Oudin", pet_Type::OUDIN},
//     {"Hamon", pet_Type::HAMON}
//   };
//   pet_Type  myPetType;
//   switch(s_mapStringToPet_Type[PetTypeStr]) {
//   case pet_Type::OUDIN:
//     myPetType = pet_Type::OUDIN;
//     // Rcpp::Rcout << "\n oudinddd \n";
//     break;
//   case pet_Type::HAMON:
//     myPetType = pet_Type::HAMON;
//     // Rcpp::Rcout << "\n hamon \n";
//     break;
//   }
//
//
//
//   hdata LatitudeVec(1,1);
//
//   LatitudeVec.resize(dHRUMdim);
//
//   for(unsigned it=0;it<dHRUMdim;it++){
//     LatitudeVec[it] = Latitude[it];
//   }
//
//   for(unsigned it=0;it<dHRUMdim;it++){
//    dHRUM_ptr.get()->calcPetToAllHrus(LatitudeVec[it], myPetType[it]);
//   }

  return ;

}

