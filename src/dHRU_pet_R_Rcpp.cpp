#include <Rcpp.h>
#include "dHRUM.h"

//' Calculates the values of Potential evapotranspiration on all singleHrus taking constant latitude
//'
//' Setting of Pet method is done by \code{PetTypeStr}, methods implemented are:
//' \code{OUDIN}, \code{HAMON}, \code{THORNTHWAITE}, \code{BLANEYCRIDDLE}, \code{JENSENHAISE}, \code{MCGUINNESSBORDNE}
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
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
//' calcPetToAllHrus(dHRUM_ptr = dhrus,50.1,"HAMON")
//' calcHBInAlldHrus(dHRUM_ptr = dhrus)
// [[Rcpp::export]]
void calcPetToAllHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, numberSel Latitude, std::string PetTypeStr){
  std::vector<std::string> PettypesVec {"OUDIN", "HAMON", "THORNTHWAITE","BLANEYCRIDDLE","JENSENHAISE", "MCGUINNESSBORDNE"};
  bool petBB=false;
  for(unsigned pp=0;pp<6;pp++){
    if(PettypesVec[pp] == PetTypeStr){
      petBB = true;
    }
  }
  if(!petBB) {
    Rcpp::stop("\n Wrong calling to the Pet Type method, try 'OUDIN' or 'HAMON or 'THORNTHWAITE' or 'BLANEYCRIDDLE' or 'JENSENHAISE' or 'MCGUINNESSBORDNE'.\n");
  }

  std::map<std::string, pet_Type> s_mapStringToPet_Type = {
    {"OUDIN", pet_Type::OUDIN},
    {"HAMON", pet_Type::HAMON},
    {"THORNTHWAITE", pet_Type::THORNTHWAITE},
    {"BLANEYCRIDDLE", pet_Type::BLANEYCRIDDLE},
    {"JENSENHAISE", pet_Type::JENSENHAISE},
    {"MCGUINNESSBORDNE", pet_Type::MCGUINNESSBORDNE}
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

//' Calculates the values of Potential evapotranspiration on all singleHrus
//'
//' Setting of Pet method is done by \code{PetTypeStr}, methods implemented are:
//' \code{OUDIN}, \code{HAMON}, \code{THORNTHWAITE}, \code{BLANEYCRIDDLE}, \code{JENSENHAISE}, \code{MCGUINNESSBORDNE}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param Latitude vector for latitude based methods
//' @param PetTypeStrNames char vector on selection of PET models
//' @param HruIds vector of Hru Ids
//' @export
//' @examples
//' nHrus <- 20
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 1,  KS = 0.01, KF = 0.03, ADIV = 0.8, CDIV = 0.3,
//'                       SDIV = 0.3, CAN_ST = 1., STEM_ST = 1., CSDIV = 0.8, TETR = 0, DDFA = 0.75, TMEL = 0.0,
//'                       RETCAP = 10 ,L = 0.5, D_BYPASS = 0.8, B_EXP=1, KS2 = 0.1, THR = 10, ALPHA =0.5)
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
//' latits = runif(nHrus)
//' pets =c("OUDIN", "HAMON", "THORNTHWAITE","BLANEYCRIDDLE","JENSENHAISE", "MCGUINNESSBORDNE")
//' PetsNams = rep(pets,times=4)
//' calcPetToAllHrusDist(dHRUM_ptr = dhrus,latits,PetTypeStrNames = PetsNams[1:nHrus],IdsHrus)
//  [[Rcpp::export]]
void calcPetToAllHrusDist(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector Latitude, Rcpp::CharacterVector PetTypeStrNames, Rcpp::CharacterVector HruIds){
  std::vector<std::string>  petNameStr = Rcpp::as<std::vector<std::string> >(PetTypeStrNames);
  unsigned dHRUMdim = dHRUM_ptr->getdHRUdim();

  for(unsigned it=0; it<dHRUMdim;it++ ){
    if ( std::find(allPetNames.begin(), allPetNames.end(), petNameStr[it]) == allPetNames.end()) {
      Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
      Rcpp::Rcout << "\n Wrong names of Pet names.\n";
      Rcpp::stop("\n Wrong calling to the Pet Type method, try:\n 'OUDIN' nor 'HAMON or 'THORNTHWAITE'\n or 'BLANEYCRIDDLE' or 'JENSENHAISE' or 'MCGUINNESSBORDNE'.\n");
    }
  }

  unsigned nHrusDF = Latitude.size();
  if(nHrusDF != dHRUMdim){
    Rcpp::stop("\n Different number of Hru's and Latitude values provided to dHRUM.\n It must equal to numeber of HRus.\n");
  }

  if(petNameStr.size() != dHRUMdim){
    Rcpp::stop("\n Different number of PET names and Hrus provided to dHRUM.\n It must equal to numeber of HRus.\n");
  }

  hdata LatitVecData(1,1);
  LatitVecData.resize(dHRUMdim);
  for(unsigned it=0; it<dHRUMdim;it++){
    LatitVecData[it] = Latitude[it];
    }

  dHRUM_ptr->calcPetToAllHrusDist( LatitVecData, petNameStr);

  return ;

}

