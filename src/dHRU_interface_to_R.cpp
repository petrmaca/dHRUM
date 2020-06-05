#include <Rcpp.h>
#include "dHRUM.h"

//' Initialization of dHRU pointer to a dHRU model
//'
//' Creates pointer instance of dHRU Model for the catchment.
//' initializes a dimension of dHRUM controled by the number of single Hru, areas of all single HRU's,
//' and ID's of all single HRU units.
//'
//' @param dimdHru a single \code{numberDta} number of single HRU units.
//' @param vecAreas a \code{numeric vector} of size \code{dimHru} of Areas for all single HRUs on dHRU.
//' @param hrusIDs a \code{character vector} of size \code{dimHru} of Id's for all single HRUs on dHRU.
//' @return dHRU_ptr pointer to dHru instance.
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
// [[Rcpp::export]]
Rcpp::XPtr<dHRUM> initdHruModel(numberDta dimdHru, Rcpp::NumericVector vecAreas, Rcpp::StringVector hrusIDs) {
  // implicit form
  // 1) creates dHRU instance on the heap (allocates memory and call constructor with no arguments)
  // 2) creates dhruPTR variable on the stack initialized with pointer to dHRU's instance
  dHRUM* dHRUM_ptr = new dHRUM();

  unsigned vecArSize = vecAreas.size();
  unsigned vecIdNamesSize = hrusIDs.size();
  if(!((dimdHru == vecArSize) && (dimdHru == vecIdNamesSize) && (vecIdNamesSize == vecIdNamesSize))) {
    Rcpp::stop("\nThe dim of dHRU single HRU units does not correspond to the length of Areas or Ids.\n");
  } else {
    single_HMunit sHRU_to_VEC;
    dHRUM_ptr->initHrusVec(dimdHru, sHRU_to_VEC);
    std::vector<std::string> vecIDs;
    vecIDs.resize(dimdHru);
    for(unsigned id=0;id<dimdHru;id++){
      vecIDs[id] =  hrusIDs[id];
    }
    hdata vecAreasHD(1,dimdHru);
    for(unsigned aa=0;aa<dimdHru;aa++){
    vecAreasHD[aa] = (numberSel) vecAreas[aa];
    // Rcpp::Rcout << vecAreasHD[aa] << "\n";
    }
    dHRUM_ptr->initHrusID(vecIDs);
// for(unsigned aa=0;aa<dimdHru;aa++){
//   Rcpp::Rcout << "Id of  single Hru on positon " << aa << " is  " << dHRU_ptr->getSingleHruId(aa) << "\n";
// }
    dHRUM_ptr->setAreasToHrus(vecAreasHD);
  }
  return Rcpp::XPtr<dHRUM>(dHRUM_ptr);
}
//' Sets similar input data obtained from file to all single HRUs for dHRU instance.
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
//' Sets the similar values of params to dHRU model for all single HRUs.
//'
//' Setting of params to dHRU.
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @param ParsVec vector of values of parameters
//' @param ParsNames a charater vector of parameter names
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' filname2 = "../../PDM/Development/PDM_dist/data/tests/inALL/BP_1960_01_01.txt"
//' setInputsToAlldHrus(dHRUM_ptr = dhrus, filname2)
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF))
// [[Rcpp::export]]
void setParamsToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector ParsVec, Rcpp::CharacterVector ParsNames) {
  unsigned numParsNames = ParsNames.size();
  unsigned numParsVals = ParsVec.size();

  if((numParsNames!=numParsVals) || (numParsNames>15) || (numParsVals>15)) {
    Rcpp::Rcout << "The number of names of params is " << numParsNames <<"\n";
    Rcpp::stop("\n and  is diferent then required number of dHRU Par Values.\n");
  } else {
     std::vector<std::string>  parNameStr = Rcpp::as<std::vector<std::string> >(ParsNames);
     std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP","KS","KF","ADIV","CDIV", \
                                   "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR", \
                                   "DDFA","TMEL","RETCAP"};
     for(unsigned it=0; it<numParsNames;it++ ){
       if ( std::find(allParNames.begin(), allParNames.end(), parNameStr[it]) == allParNames.end()) {
        Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of Par Values.\n");
        }
      }
     std::vector<numberSel> parsVals = Rcpp::as<std::vector<numberSel>>(ParsVec);
     std::map<std::string, par_HRUtype> s_mapStringTopar_HRUtype = {
      {"B_SOIL", par_HRUtype::B_SOIL},
      {"C_MAX", par_HRUtype::C_MAX },
      {"B_EVAP",  par_HRUtype::B_EVAP},
      {"KS",  par_HRUtype::KS},
      {"KF",  par_HRUtype::KF},
      {"ADIV",  par_HRUtype::ADIV},
      {"CDIV",  par_HRUtype::CDIV},
      {"SDIV", par_HRUtype::SDIV},
      {"CAN_ST", par_HRUtype::CAN_ST},
      {"STEM_ST", par_HRUtype::STEM_ST},
      {"CSDIV", par_HRUtype::CSDIV},
      {"TETR", par_HRUtype::TETR},
      {"DDFA", par_HRUtype::DDFA},
      {"TMEL", par_HRUtype::TMEL},
      {"RETCAP", par_HRUtype::RETCAP}
      };

     std::vector<std::pair<numberSel,par_HRUtype>> ParsToLoad;

     // Rcpp::Rcout << "\n iusa " << ParsToLoad.size() << " \n";

     for(unsigned id=0;id<numParsNames;id++){
        switch(s_mapStringTopar_HRUtype[parNameStr[id]]) {
            case par_HRUtype::B_SOIL:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::B_SOIL));
              // pars[0] = par_dta;
              // pars[3] = pars[1] / (pars[0]+1);
                 // std::cout << "New b_soil --> " << parsVals[id]  << " loaded\n";
              break;
            case par_HRUtype::C_MAX:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::C_MAX));
              // TTTT pars[1] = par_dta;
              // pars[3] = pars[1] / (pars[0]+1);
                 // std::cout << "New c_max --> loaded\n";
              break;
            case par_HRUtype::B_EVAP:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::B_EVAP));
              // pars[2] = par_dta;
                 // std::cout << "New b_evap --> loaded\n";
              break;
            case par_HRUtype::SMAX:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::SMAX));
              // pars[3] = par_dta;
                 // std::cout << "New Smax --> loaded\n";
              break;
            case par_HRUtype::KS:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KS));
              // pars[4] = par_dta;
                 // std::cout << "New Ks --> loaded\n";
              break;
            case par_HRUtype::KF:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KF));
              // pars[5] = par_dta;
                 // std::cout << "New Kf --> loaded\n";
              break;
            case par_HRUtype::ADIV:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::ADIV));
              // pars[6] = par_dta;
                 // std::cout << "New Adiv --> loaded\n";
              break;
            case par_HRUtype::CDIV:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CDIV));
              // pars[7] = par_dta;
                 // std::cout << "New Cdiv --> loaded\n";
              break;
            case par_HRUtype::SDIV:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::SDIV));
              // pars[8] = par_dta;
                 // std::cout << "New Sdiv --> loaded\n";
              break;
            case par_HRUtype::CAN_ST:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CAN_ST));
              // pars[9] = par_dta;
                 // std::cout << "New Can_St --> loaded\n";
              break;
            case par_HRUtype::STEM_ST:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::STEM_ST));
              // pars[10] = par_dta;
                 // std::cout << "New Stem_St --> loaded\n";
              break;
            case par_HRUtype::CSDIV:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CSDIV));
              // pars[11] = par_dta;
                 // std::cout << "New CSdiv --> loaded\n";
              break;
            case par_HRUtype::TETR:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::TETR));
              // pars[12] = par_dta;
                 // std::cout << "New TETR --> loaded\n";
              break;
            case par_HRUtype::DDFA:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::DDFA));
              // pars[13] = par_dta;
                 // std::cout << "New DDFA --> loaded\n";
              break;
            case par_HRUtype::TMEL:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::TMEL));
              // pars[14] = par_dta;
                 // std::cout << "New TMEL --> loaded\n";
              break;
            case par_HRUtype::RETCAP:
              ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RETCAP));
              // pars[15] = par_dta;
                 // std::cout << "New RETCAP --> loaded\n";
              break;
            }
        }
     dHRUM_ptr.get()->setParamsToAllHrus(ParsToLoad);
    }
 return  ;
}
//' Calculates the values of Potetial evpotranspiration on all singleHrus
//'
//' Setting of Pet method is done by \code{PetTypeStr}, methods implemented are: \code{Oudin}, \code{Hamon}
//'
//'
//' @param dHRU_ptr pointer to dHRU instance
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
//' Calculates the values of fluxes and state variables for all Hrus in dHru
//'
//' Updating the states ad values of all ts varibles in all single Hrus of dHRU
//'
//'
//' @param dHRU_ptr pointer to dHRU instance
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
//' filname2 = "../data/inBP_1960_01_01.txt"
//' setInputsToAlldHrus(dhrus, filname2)
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dhrus,as.numeric(ParDF[1,]),names(ParDF),TRUE,0)
//' calcPetToHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
//' calcHBInAlldHrus(dHRUM_ptr = dhrus)
//' gatherHBdata(dHRUM_ptr = dhrus)
//' printToFile(dhrus,file)
// [[Rcpp::export]]
void printToFile(Rcpp::XPtr<dHRUM> dHRUM_ptr, std::string namOutFilet) {

  return dHRUM_ptr.get()->printAllDta(namOutFilet);

}
//' Provides with dHRU outputs on ts
//'
//' return matrix with state varibles and fluxes averadged over basin area.
//'
//'
//' @param dHRU_ptr pointer to dHRU instance
//' @return a list with matrix of caldata \code{[,1:4]} and hdata ts variables \code{[,5:27]} and names of vars
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' filname2 = "../tests/indata/BP_1960_01_01.txt"
//' setInputsToAlldHrus(dHRUM_ptr = dhrus, filname2)
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF),TRUE,0)
//' calcPetToHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
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
  Rcpp::StringVector VarsNams({"YEAR", "MONTH", "DAY", "JDAY", \
   "PREC","SNOW","AET","PET","TEMP", \
   "MELT","TROF","STEF","CANF","CANS", \
   "STES","EVAC","EVAS","EVBS","INTS", \
   "SOIS","GROS","SURS","TOTR","BASF","DIRR","PERC","PREF"});

  return Rcpp::List::create(
    Rcpp::Named("outDta") = outDta,
    Rcpp::Named("VarsNams") = VarsNams
    );
}

//' Sets the values of params to dHRU model for one single HRU.
//'
//' Setting of params to one particular single HRU for dHRU.
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @param ParsVec vector of values of parameters
//' @param ParsNames a charater vector of parameter names
//' @param singleHruId the numerical single HRU ID
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
//' setParamsToOnedHru(dHRUM_ptr = dhrus,as.numeric(ParDF[1,]),names(ParDF),0)
// [[Rcpp::export]]
void setParamsToOnedHru(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector ParsVec, Rcpp::CharacterVector ParsNames, unsigned singleHruId) {
  unsigned numParsNames = ParsNames.size();
  unsigned numParsVals = ParsVec.size();

  if((numParsNames!=numParsVals) || (numParsNames>15) || (numParsVals>15)) {
    Rcpp::Rcout << "The number of names of params is " << numParsNames <<"\n";
    Rcpp::stop("\n and  is diferent then required number of dHRU Par Values.\n");
  } else {
    std::vector<std::string>  parNameStr = Rcpp::as<std::vector<std::string> >(ParsNames);
    std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP","KS","KF","ADIV","CDIV", \
                                          "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR", \
                                          "DDFA","TMEL","RETCAP"};
    for(unsigned it=0; it<numParsNames;it++ ){
      if ( std::find(allParNames.begin(), allParNames.end(), parNameStr[it]) == allParNames.end()) {
        Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of Par Values.\n");
      }
    }
    std::vector<numberSel> parsVals = Rcpp::as<std::vector<numberSel>>(ParsVec);
    std::map<std::string, par_HRUtype> s_mapStringTopar_HRUtype = {
      {"B_SOIL", par_HRUtype::B_SOIL},
      {"C_MAX", par_HRUtype::C_MAX },
      {"B_EVAP",  par_HRUtype::B_EVAP},
      {"KS",  par_HRUtype::KS},
      {"KF",  par_HRUtype::KF},
      {"ADIV",  par_HRUtype::ADIV},
      {"CDIV",  par_HRUtype::CDIV},
      {"SDIV", par_HRUtype::SDIV},
      {"CAN_ST", par_HRUtype::CAN_ST},
      {"STEM_ST", par_HRUtype::STEM_ST},
      {"CSDIV", par_HRUtype::CSDIV},
      {"TETR", par_HRUtype::TETR},
      {"DDFA", par_HRUtype::DDFA},
      {"TMEL", par_HRUtype::TMEL},
      {"RETCAP", par_HRUtype::RETCAP}
    };

    std::vector<std::pair<numberSel,par_HRUtype>> ParsToLoad;

    // Rcpp::Rcout << "\n iusa " << ParsToLoad.size() << " \n";

    for(unsigned id=0;id<numParsNames;id++){
      switch(s_mapStringTopar_HRUtype[parNameStr[id]]) {
      case par_HRUtype::B_SOIL:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::B_SOIL));
        // pars[0] = par_dta;
        // pars[3] = pars[1] / (pars[0]+1);
        // std::cout << "New b_soil --> " << parsVals[id]  << " loaded\n";
        break;
      case par_HRUtype::C_MAX:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::C_MAX));
        // TTTT pars[1] = par_dta;
        // pars[3] = pars[1] / (pars[0]+1);
        // std::cout << "New c_max --> loaded\n";
        break;
      case par_HRUtype::B_EVAP:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::B_EVAP));
        // pars[2] = par_dta;
        // std::cout << "New b_evap --> loaded\n";
        break;
      case par_HRUtype::SMAX:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::SMAX));
        // pars[3] = par_dta;
        // std::cout << "New Smax --> loaded\n";
        break;
      case par_HRUtype::KS:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KS));
        // pars[4] = par_dta;
        // std::cout << "New Ks --> loaded\n";
        break;
      case par_HRUtype::KF:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KF));
        // pars[5] = par_dta;
        // std::cout << "New Kf --> loaded\n";
        break;
      case par_HRUtype::ADIV:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::ADIV));
        // pars[6] = par_dta;
        // std::cout << "New Adiv --> loaded\n";
        break;
      case par_HRUtype::CDIV:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CDIV));
        // pars[7] = par_dta;
        // std::cout << "New Cdiv --> loaded\n";
        break;
      case par_HRUtype::SDIV:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::SDIV));
        // pars[8] = par_dta;
        // std::cout << "New Sdiv --> loaded\n";
        break;
      case par_HRUtype::CAN_ST:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CAN_ST));
        // pars[9] = par_dta;
        // std::cout << "New Can_St --> loaded\n";
        break;
      case par_HRUtype::STEM_ST:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::STEM_ST));
        // pars[10] = par_dta;
        // std::cout << "New Stem_St --> loaded\n";
        break;
      case par_HRUtype::CSDIV:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CSDIV));
        // pars[11] = par_dta;
        // std::cout << "New CSdiv --> loaded\n";
        break;
      case par_HRUtype::TETR:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::TETR));
        // pars[12] = par_dta;
        // std::cout << "New TETR --> loaded\n";
        break;
      case par_HRUtype::DDFA:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::DDFA));
        // pars[13] = par_dta;
        // std::cout << "New DDFA --> loaded\n";
        break;
      case par_HRUtype::TMEL:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::TMEL));
        // pars[14] = par_dta;
        // std::cout << "New TMEL --> loaded\n";
        break;
      case par_HRUtype::RETCAP:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RETCAP));
        // pars[15] = par_dta;
        // std::cout << "New RETCAP --> loaded\n";
        break;
      }
    }

    unsigned dHRUdim = 0;
    dHRUdim = dHRUM_ptr.get()->getdHRUdim();
    if( singleHruId > dHRUdim){
      Rcpp::stop("The wrong ID or number of single HRU.\n");
    }

    dHRUM_ptr.get()->setParamsToOneHru(ParsToLoad,singleHruId);
  }
  return  ;
}
//' Sets the Precipitation and temperature vectors to dHRU.
//'
//' Setting the similar vecotr of Precipitation and temperature to all single HRU.
//' Setting the calender using the first date fo period
//'
//' @param dHRUM_ptr pointer to dHRU instance
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
