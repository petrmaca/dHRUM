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
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(ParDF[1,]),ParsNames =names(ParDF))
// [[Rcpp::export]]
void setParamsToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::NumericVector ParsVec, Rcpp::CharacterVector ParsNames) {
  unsigned numParsNames = ParsNames.size();
  unsigned numParsVals = ParsVec.size();

  if((numParsNames!=numParsVals) || (numParsNames>33) || (numParsVals>33)) {
    Rcpp::Rcout << "The number of names of params is " << numParsNames <<"\n";
    Rcpp::Rcout << "The number of values of params is " << numParsVals <<"\n";
    Rcpp::stop("\n Those values are different or higher or smaller then required number of dHRU Par Values =22.\n");
  } else {

    std::vector<std::string>  parNameStr = Rcpp::as<std::vector<std::string> >(ParsNames);

    std::cout << "Pars names prints \n";
    for(unsigned i=0;i<numParsNames;i++){
      std::cout << parNameStr[i] << std::endl;
    }
    std::cout << "\n";
/*
    // std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP","KS","KF","ADIV","CDIV", \
    //                               "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR",         \
    //                               "DDFA","TMEL","RETCAP"};

 */
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
      {"SMAX", par_HRUtype::SMAX },
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
      {"RETCAP", par_HRUtype::RETCAP},
      {"L", par_HRUtype::L},
      {"D_BYPASS", par_HRUtype::D_BYPASS},
      {"B_EXP", par_HRUtype::B_EXP},
      {"KS2", par_HRUtype::KS2},
      {"THR", par_HRUtype::THR},
      {"ALPHA", par_HRUtype::ALPHA},
      {"CMIN", par_HRUtype::CMIN},
      {"FC", par_HRUtype::FC},
      {"FOREST_FRACT", par_HRUtype::FOREST_FRACT},
      {"KF2", par_HRUtype::KF2},
      {"KF_NONLIN", par_HRUtype::KF_NONLIN},
      {"C", par_HRUtype::C},
      {"INFR_MAX", par_HRUtype::INFR_MAX},
      {"RF", par_HRUtype::RF},
      {"WP", par_HRUtype::WP},
      {"RBAI", par_HRUtype::RBAI},
      {"RBEI", par_HRUtype::RBEI},
      {"KFR", par_HRUtype::KFR}

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
      case par_HRUtype::SMAXpdm:

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
      case par_HRUtype::L:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::L));
        // pars[16] = par_dta;
        // std::cout << "New L --> loaded\n";
        break;
      case par_HRUtype::D_BYPASS:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::D_BYPASS));
        // pars[17] = par_dta;
        // std::cout << "New D_BYPASS --> loaded\n";
        break;
      case par_HRUtype::B_EXP:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::B_EXP));
        // pars[18] = par_dta;
        // std::cout << "New B_EXP --> loaded\n";
        break;
      case par_HRUtype::KS2:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KS2));
        // pars[19] = par_dta;
        // std::cout << "New KS2 --> loaded\n";
        break;
      case par_HRUtype::THR:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::THR));
        // pars[20] = par_dta;
        // std::cout << "New THR --> loaded\n";
        break;
      case par_HRUtype::ALPHA:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::ALPHA));
        // pars[21] = par_dta;
        // std::cout << "New ALPHA --> loaded\n";
        break;
      case par_HRUtype::CMIN:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CMIN));
        // pars[22] = par_dta;
        // std::cout << "New ALPHA --> loaded\n";
        break;
      case par_HRUtype::FC:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::FC));
        // pars[23] = par_dta;
        // std::cout << "New FC --> loaded\n";
        break;
      case par_HRUtype::FOREST_FRACT:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::FOREST_FRACT));
        // pars[24] = par_dta;
        // std::cout << "New FOREST_FRACT --> loaded\n";
        break;
      case par_HRUtype::KF2:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KF2));
        // pars[25] = par_dta;
        // std::cout << "New KF2 --> loaded\n";
        break;
      case par_HRUtype::KF_NONLIN:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KF_NONLIN));
        // pars[26] = par_dta;
        // std::cout << "New KF_NONLIN --> loaded\n";
        break;
      case par_HRUtype::C:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::C));
        // pars[27] = par_dta;
        // std::cout << "New C --> loaded\n";
        break;
      case par_HRUtype::INFR_MAX:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::INFR_MAX));
        // pars[28] = par_dta;
        // std::cout << "New INFR_MAX --> loaded\n";
        break;
      case par_HRUtype::RF:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RF));
        // pars[29] = par_dta;
        // std::cout << "New RF --> loaded\n";
        break;
      case par_HRUtype::WP:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::WP));
        // pars[30] = par_dta;
        // std::cout << "New WP --> loaded\n";
      case par_HRUtype::RBAI:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RBAI));
        // pars[31] = par_dta;
        // std::cout << "New RBAI --> loaded\n";
        break;
      case par_HRUtype::RBEI:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RBEI));
        // pars[32] = par_dta;
        // std::cout << "New RBEI --> loaded\n";
      case par_HRUtype::KFR:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KFR));
        // pars[34] = par_dta;
        // std::cout << "New KFR --> loaded\n";

        break;
      }
    }

    ///neuplna funkce getRequiredParamsForHru(id)?????
    unsigned dimHRU = dHRUM_ptr.get()->getdHRUdim();
    for(unsigned id=0;id<dimHRU;id++){
      std::vector<std::string> requiredParams = dHRUM_ptr.get()->getRequiredParamsForHru(id);
      for(unsigned i=0;i<requiredParams.size();i++){
        if ( std::find(ParsNames.begin(), ParsNames.end(), requiredParams[i]) == ParsNames.end()) {
          Rcpp::Rcout << "\nSomething wrong in Pars setting on item " << (i+1) << "\n";
          Rcpp::stop("\n Required parameter is missing.\n");
        }
      }
    }
    dHRUM_ptr.get()->setParamsToAlldHrus(ParsToLoad);
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
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
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

/*
//   std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP","KS","KF","ADIV","CDIV", \
//                                          "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR", \
//                                          "DDFA","TMEL","RETCAP"};
*/
     for(unsigned it=0; it<numParsNames;it++ ){
      if ( std::find(allParNames.begin(), allParNames.end(), parNameStr[it]) == allParNames.end()) {
        Rcpp::Rcout << "\nSomething wrong on for Params item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of Par Values.\n");
      }
    }
    std::vector<numberSel> parsVals = Rcpp::as<std::vector<numberSel>>(ParsVec);
    std::map<std::string, par_HRUtype> s_mapStringTopar_HRUtype = {
      {"B_SOIL", par_HRUtype::B_SOIL},
      {"C_MAX", par_HRUtype::C_MAX },
      {"SMAX", par_HRUtype::SMAX },
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
      {"RETCAP", par_HRUtype::RETCAP},
      {"L", par_HRUtype::L},
      {"D_BYPASS", par_HRUtype::D_BYPASS},
      {"B_EXP", par_HRUtype::B_EXP},
      {"KS2", par_HRUtype::KS2},
      {"THR", par_HRUtype::THR},
      {"ALPHA", par_HRUtype::ALPHA},
      {"FC", par_HRUtype::FC},
      {"FOREST_FRACT", par_HRUtype::FOREST_FRACT},
      {"KF2", par_HRUtype::KF2},
      {"KF_NONLIN", par_HRUtype::KF_NONLIN},
      {"C", par_HRUtype::C},
      {"INFR_MAX", par_HRUtype::INFR_MAX},
      {"RF", par_HRUtype::RF},
      {"WP", par_HRUtype::WP}
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
      case par_HRUtype::SMAXpdm:

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
      case par_HRUtype::L:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::L));
        // pars[16] = par_dta;
        // std::cout << "New L --> loaded\n";
        break;
      case par_HRUtype::D_BYPASS:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::D_BYPASS));
        // pars[17] = par_dta;
        // std::cout << "New D_BYPASS --> loaded\n";
        break;
      case par_HRUtype::B_EXP:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::B_EXP));
        // pars[18] = par_dta;
        // std::cout << "New B_EXP --> loaded\n";
        break;
      case par_HRUtype::KS2:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KS2));
        // pars[19] = par_dta;
        // std::cout << "New KS2 --> loaded\n";
        break;
      case par_HRUtype::THR:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::THR));
        // pars[20] = par_dta;
        // std::cout << "New THR --> loaded\n";
        break;
      case par_HRUtype::ALPHA:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::ALPHA));
        // pars[21] = par_dta;
        // std::cout << "New ALPHA --> loaded\n";
        break;
      case par_HRUtype::CMIN:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::CMIN));
        // pars[22] = par_dta;
        // std::cout << "New ALPHA --> loaded\n";
        break;
      case par_HRUtype::FC:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::FC));
        // pars[23] = par_dta;
        // std::cout << "New FC --> loaded\n";
        break;
      case par_HRUtype::FOREST_FRACT:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::FOREST_FRACT));
        // pars[24] = par_dta;
        // std::cout << "New FOREST_FRACT --> loaded\n";
        break;
      case par_HRUtype::KF2:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KF2));
        // pars[25] = par_dta;
        // std::cout << "New KF2 --> loaded\n";
        break;
      case par_HRUtype::KF_NONLIN:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KF_NONLIN));
        // pars[26] = par_dta;
        // std::cout << "New KF_NONLIN --> loaded\n";
        break;
      case par_HRUtype::C:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::C));
        // pars[27] = par_dta;
        // std::cout << "New C --> loaded\n";
        break;
      case par_HRUtype::INFR_MAX:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::INFR_MAX));
        // pars[28] = par_dta;
        // std::cout << "New INFR_MAX --> loaded\n";
        break;
      case par_HRUtype::RF:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RF));
        // pars[29] = par_dta;
        // std::cout << "New RF --> loaded\n";
        break;
      case par_HRUtype::WP:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::WP));
        // pars[30] = par_dta;
        // std::cout << "New WP --> loaded\n";
        break;
      case par_HRUtype::RBAI:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RBAI));
        // pars[31] = par_dta;
        // std::cout << "New RBAI --> loaded\n";
        break;
      case par_HRUtype::RBEI:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::RBEI));
        // pars[32] = par_dta;
        // std::cout << "New RBEI --> loaded\n";
        break;
      case par_HRUtype::KFR:
        ParsToLoad.push_back(std::make_pair((numberSel) parsVals[id], par_HRUtype::KFR));
        // pars[32] = par_dta;
        // std::cout << "New RBEI --> loaded\n";
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

//' Sets the values of params to distributed dHRU model using data.frame of params as inputs.
//'
//' Setting of vectors of params to all HRUs for distributed dHRUM.
//'
//' @param dHRUM_ptr pointer to dHRU instance
//' @param ParsDF data.frame of parametrs cols show parameters, rows show the Hrus
//' @param PrintPars if TRUE than params are printed
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//' DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' parsvec=as.numeric(ParDF[1,])
//' for(i in 1:(nHrus-1)){
//' ParDF <- rbind(ParDF,parsvec)
//' }
//' ParDF[,1] <- seq(1:nHrus)
//' setParsToDistdHRUM(dhrus, ParDF, FALSE)
//' ParDF[,2] <- seq(1:nHrus)
//' setParsToDistdHRUM(dhrus, ParDF, FALSE)
// [[Rcpp::export]]
void setParsToDistdHRUM(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::DataFrame ParsDF, bool PrintPars) {

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
      Rcpp::Rcout << "\nSomething wrong on parameter names in item " << (it+1) << "\n";
      Rcpp::stop("\nWrong names of Par Values.\n");
    }
  }
// #pragma omp parallel for
  for(unsigned it=0;it<dimDHRUM; it++){
    // Rcpp::Rcout << it << "\n";
    Rcpp::NumericVector ParsVec(numColsParsMat);
    for(unsigned cl=0;cl<numColsParsMat;cl++){
     Rcpp::NumericVector ParVectsCols =  ParsDF[cl];
      ParsVec[cl] = ParVectsCols[it];
    }
    setParamsToOnedHru( dHRUM_ptr, ParsVec, ParsNames, it);
  }

  if(PrintPars) dHRUM_ptr.get()->print_Pars();

  return  ;
}


//' Getting the current singeHMunit parameters.
//'
//' shows the DataFrame of parameters for selected HRU
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param singleHruId a Id of particular Hru
//' @export
//' @examples
//' nHrus <- 1
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(ParDF[1,]),ParsNames =names(ParDF))
//' getCurdHRUpars(dHRUM_ptr = dhrus,0)
// [[Rcpp::export]]
 Rcpp::DataFrame getCurdHRUpars(Rcpp::XPtr<dHRUM> dHRUM_ptr,unsigned singleHruId) {


   Rcpp::NumericVector cur_par;
   Rcpp::NumericVector up_par;
   Rcpp::NumericVector low_par;
   Rcpp::StringVector par_names;
   dHRUM_ptr.get()->Current_Params(singleHruId);
   cur_par = dHRUM_ptr.get()->get_param_vec(singleHruId);
   up_par = dHRUM_ptr.get()->get_upparam_vec(singleHruId);
   low_par = dHRUM_ptr.get()->get_lowparam_vec(singleHruId);
   par_names = dHRUM_ptr.get()->get_param_names(singleHruId);

   return Rcpp::DataFrame::create(
     Rcpp::Named("Cur_names") = par_names,
     Rcpp::Named("Cur_par") = cur_par,
     Rcpp::Named("Up_bound") = up_par,
     Rcpp::Named("Low_bound") = low_par);
 }

//' Getting the all HM units parameters.
 //'
 //' shows the list of parameters for all HRUs
 //'
 //' @param dHRUM_ptr pointer to dHRUM instance
 //' @export
 //' @examples
 //' nHrus <- 1
 //' Areas <- runif(nHrus,min = 1,max  = 10)
 //' IdsHrus <- paste0("ID",seq(1:length(Areas)))
 //' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
 //' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
 //' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
 //' prec=c(1,2,3)
 //' temp=c(1,2,3)
 //' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
 //'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
 //' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
 //'  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
 //' setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(ParDF[1,]),ParsNames =names(ParDF))
 //' getAllHRUpars(dHRUM_ptr = dhrus)
 // [[Rcpp::export]]
 Rcpp::List getAllHRUpars(Rcpp::XPtr<dHRUM> dHRUM_ptr) {

   Rcpp::DataFrame df;
   Rcpp::List alldfs;

   int HRUnum = dHRUM_ptr.get()-> getdHRUdim(); //pocet hru
   Rcpp::NumericVector cur_par;
   Rcpp::NumericVector up_par;
   Rcpp::NumericVector low_par;
   Rcpp::StringVector par_names;
   Rcpp::StringVector HRU_name;

   for ( int HruId = 0; HruId< HRUnum; HruId++ )
   {
     dHRUM_ptr.get()->Current_Params(HruId);
     cur_par = dHRUM_ptr.get()->get_param_vec(HruId);
     up_par = dHRUM_ptr.get()->get_upparam_vec(HruId);
     low_par = dHRUM_ptr.get()->get_lowparam_vec(HruId);
     par_names = dHRUM_ptr.get()->get_param_names(HruId);
     std::string textID = dHRUM_ptr.get()->getSingleHruId(HruId);
     HRU_name=textID;

     df=Rcpp::DataFrame::create(
       Rcpp::Named("Cur_names") = par_names,
       Rcpp::Named("Cur_par") = cur_par,
       Rcpp::Named("Up_bound") = up_par,
       Rcpp::Named("Low_bound") = low_par,
       Rcpp::Named("sHRU_ID") = HRU_name);

     alldfs.insert( HruId, df );
   }

   return alldfs;
 }

//' Getting the current singeHMunit configuration.
//'
//' shows the list of configuration for selected HRU
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param singleHruId a Id of particular Hru
//' @export
//' @examples
//' nHrus <- 1
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(ParDF[1,]),ParsNames =names(ParDF))
//' getCurSHRUconfig(dHRUM_ptr = dhrus,0)
// [[Rcpp::export]]
Rcpp::DataFrame getCurSHRUconfig(Rcpp::XPtr<dHRUM> dHRUM_ptr,unsigned singleHruId) {

   Rcpp::StringVector names;
   Rcpp::StringVector values;
   Rcpp::DataFrame df;

   auto HRUnum = dHRUM_ptr.get()-> getdHRUdim(); //pocet hru

   if(singleHruId>(HRUnum-1)){
     std::cout<<"Wrong sHRU value!! Currently exist "<< HRUnum<<" sHRUs! Indexing starts from 0" <<std::endl;
   }else {
     std::vector<std::pair<std::string,std::string>> sHruConfig;
     sHruConfig=dHRUM_ptr.get()->get_sHMu_Config(singleHruId);

     for ( auto it = sHruConfig.begin(); it != sHruConfig.end(); it++ )
     {
       names.push_back(it->first);
       values.push_back(it->second);
     }

     std::string textID = dHRUM_ptr.get()->getSingleHruId(singleHruId);
     df = Rcpp::DataFrame::create( Rcpp::Named("V1") = names,
                                Rcpp::Named(textID) = values);
   }

   return df;
}


//' Getting configurations of all HRUs.
//'
//' shows the data frame of all HRUs configurations
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @export
//' @examples
//' nHrus <- 10
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
//' setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' prec=c(1,2,3)
//' temp=c(1,2,3)
//' setPTDateInputsToAlldHrus(dhrus, Prec = prec, Temp = temp,
//'   DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//'  SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' setParamsToAlldHrus(dHRUM_ptr = dhrus,ParsVec = as.numeric(ParDF[1,]),ParsNames =names(ParDF))
//' getAllHRUconfigs(dHRUM_ptr = dhrus,IdsHrus)
// [[Rcpp::export]]
Rcpp::DataFrame getAllHRUconfigs(Rcpp::XPtr<dHRUM> dHRUM_ptr) {

   Rcpp::StringVector names;
   std::vector<std::pair<std::string,std::string>> sHruConfig;

   int HRUnum = dHRUM_ptr.get()-> getdHRUdim(); //pocet hru

   sHruConfig=dHRUM_ptr.get()->get_sHMu_Config(0);
   for ( auto it = sHruConfig.begin(); it != sHruConfig.end(); it++ )
   {     names.push_back(it->first);  }

   Rcpp::DataFrame df=Rcpp::DataFrame::create( Rcpp::Named("Structure") = names);

   for ( int i = 0; i< HRUnum; i++ )
   {
     sHruConfig=dHRUM_ptr.get()->get_sHMu_Config(i);
     Rcpp::StringVector values;
     for ( auto it = sHruConfig.begin(); it != sHruConfig.end(); it++ )
     {
       values.push_back(it->second);
     }
     std::string textID = dHRUM_ptr.get()->getSingleHruId(i);
     df.push_back(values,textID);
    }

   return df;
 }
