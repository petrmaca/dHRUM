#include <Rcpp.h>
#include "dHRUM.h"

//' Sets the similar values of params to dHRU model for all single HRUs.
//'
//' Setting of params to dHRUM.
//'
//' @param dHRUM_ptr pointer to dHRUM instance
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
    // std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP","KS","KF","ADIV","CDIV", \
    //                               "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR",         \
    //                               "DDFA","TMEL","RETCAP"};
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
  unsigned numPars1Hru;

  numPars1Hru = dHRUM_ptr.get()->get_singleHRUnumPars(singleHruId);

  if((numParsNames!=numParsVals) || (numParsNames>numPars1Hru) || (numParsVals>numPars1Hru)) {
    Rcpp::Rcout << "The number of names of params is " << numParsNames <<"\n";
    Rcpp::stop("\n and  is diferent then required number of dHRU Par Values.\n");
  } else {
    std::vector<std::string>  parNameStr = Rcpp::as<std::vector<std::string> >(ParsNames);
//   std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP","KS","KF","ADIV","CDIV", \
//                                          "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR", \
//                                          "DDFA","TMEL","RETCAP"};
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

//' Sets the values of params to dHRU model using data.frame of params as inputs.
//'
//' Setting of vectors of params to all HRUs for distributed dHRUM.
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @param ParsDF data.frame of parametrs cols show parameters, rows show the Hrus
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
void setParsToDistdHRUM(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::DataFrame ParsDF) {

  unsigned dimDHRUM=0, numColsParsMat=0;
  dimDHRUM = dHRUM_ptr.get()->getdHRUdim();//Should it be added into the dHRUM class??
  numColsParsMat = ParsDF.length();
  // Rcpp::Rcout <<dimDHRUM << " dimhru " << numColsParsMat << "nomclos \n";
  //ToDo check on number of cols for ParsDF smaller and equal to number of pars
  // and bigger than 0
  // checking the correct number of params to load
  for(unsigned it=0;it<dimDHRUM; it++){
    if(numColsParsMat > dHRUM_ptr.get()->get_singleHRUnumPars(it)){
      Rcpp::Rcout << "The single HRU with Id " << dHRUM_ptr.get()->getSingleHruId(it) << " have different number of params.\n";
      Rcpp::stop("\nWrong size number of parameters.\n");
    }
  }

  Rcpp::CharacterVector ParsNames = ParsDF.names();
  ParsNames = ParsDF.names();
  std::vector<std::string>  parNameStr = Rcpp::as<std::vector<std::string> >(ParsNames);
  for(unsigned it=0; it<numColsParsMat;it++ ){
    if ( std::find(allParNames.begin(), allParNames.end(), parNameStr[it]) == allParNames.end()) {
      Rcpp::Rcout << "\nSomething wrong parameter names on item " << (it+1) << "\n";
      Rcpp::stop("\nWrong names of Par Values.\n");
    }
  }

  for(unsigned it=0;it<dimDHRUM; it++){
    Rcpp::NumericVector ParsVec(numColsParsMat);
    for(unsigned cl=0;cl<numColsParsMat;cl++){
     Rcpp::NumericVector ParVectsCols =  ParsDF[cl];
      ParsVec[cl] = ParVectsCols[it];
    }
    setParamsToOnedHru( dHRUM_ptr, ParsVec, ParsNames, it);
  }

  return  ;
}
