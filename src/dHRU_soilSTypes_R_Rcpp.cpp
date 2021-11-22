#include <Rcpp.h>
#include "dHRUM.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
void setSoilStorTypeToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector soilTypes, Rcpp::CharacterVector hruIds) {

  unsigned numSoilSTypes = soilTypes.size();
  unsigned numHruIdNames = hruIds.size();

  //check if names are consistent
  //for which hrus we want to change the stor type - vector of character
  if(numSoilSTypes!=numHruIdNames) {
    Rcpp::Rcout << "The number of soil types does not correspond to the number of HRUs " << numSoilSTypes <<"\n";
    Rcpp::stop("\nWrong size number of soilS types.\n");
  } else {
    std::vector<std::string> soilNameStr = Rcpp::as<std::vector<std::string> >(soilTypes);
    for(unsigned it=0; it<numSoilSTypes;it++ ){
      if ( std::find(allSoilStorTypeNames.begin(), allSoilStorTypeNames.end(), soilNameStr[it]) == allSoilStorTypeNames.end()) {
        Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of soilS type values.\n");
      }
    }
    const std::vector<std::string> ids = dHRUM_ptr.get()->getHRUIds();
    std::vector<std::string> hruIdName = Rcpp::as<std::vector<std::string> >(hruIds);

    //std::vector<std::string> ids = Rcpp::as<std::vector<std::string> >(hruIds);
    for(unsigned it=0; it<numHruIdNames;it++ ){
      if ( std::find(ids.begin(), ids.end(), hruIdName[it]) == ids.end()) {
        Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of Hru Id Values.\n");
      }
    }
    std::map<std::string, soil_STORtype> s_mapStringToSoilSType_HRUtype = {
      {"PDM", soil_STORtype::PDM},
      {"COLLIE_V2", soil_STORtype::COLLIE_V2},
      {"NEW_ZEALAND", soil_STORtype::NEW_ZEALAND},
      {"GR4J", soil_STORtype::GR4J},
      {"SBROOK_V1", soil_STORtype::SBROOK_V1},
      {"HILLSLOPE", soil_STORtype::HILLSLOPE},
      {"PLATEAU", soil_STORtype::PLATEAU}
    };
    std::vector<unsigned> indexHru;
    indexHru.resize(hruIds.size());
    for(unsigned i=0; i<indexHru.size(); i++) {
      for(unsigned j=0; j<dHRUM_ptr.get()->getdHRUdim(); j++) {
        if(!hruIdName[i].compare(ids[j])) {
          indexHru[i] = j;
        }
      }
      // std::cout << indexHru[i];
    }
    std::vector<std::pair<unsigned,soil_STORtype>> soilTypesToLoad;


    for(unsigned id=0;id<numHruIdNames;id++) {
      switch(s_mapStringToSoilSType_HRUtype[soilNameStr[id]]) {
      case soil_STORtype::PDM:
        soilTypesToLoad.push_back(std::make_pair(indexHru[id], soil_STORtype::PDM));
        break;
      case soil_STORtype::COLLIE_V2:
        soilTypesToLoad.push_back(std::make_pair(indexHru[id], soil_STORtype::COLLIE_V2));
        break;
      case soil_STORtype::NEW_ZEALAND:
        soilTypesToLoad.push_back(std::make_pair(indexHru[id], soil_STORtype::NEW_ZEALAND));
        break;
      case soil_STORtype::GR4J:
        soilTypesToLoad.push_back(std::make_pair(indexHru[id], soil_STORtype::GR4J));
        break;
      case soil_STORtype::SBROOK_V1:
        soilTypesToLoad.push_back(std::make_pair(indexHru[id], soil_STORtype::SBROOK_V1));
        break;
      case soil_STORtype::HILLSLOPE:
        soilTypesToLoad.push_back(std::make_pair(indexHru[id], soil_STORtype::HILLSLOPE));
        break;
      case soil_STORtype::PLATEAU:
        soilTypesToLoad.push_back(std::make_pair(indexHru[id], soil_STORtype::PLATEAU));
        break;
      }
    }
    dHRUM_ptr.get()->initSoilStypeToAlldHrus(soilTypesToLoad);
  }

  return;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
