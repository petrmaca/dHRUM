#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of interception models types to dHRU model for all single HRUs.
//'
//' Setting the interception type to dHRUM to all HRUs. Possibe types: \code{Rutter_Gash}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param intcptnTypes a charater vector of Interception type names
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
void setInterceptiontypeToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector intcptnTypes, Rcpp::CharacterVector hruIds) {
  unsigned numINTRTypes = intcptnTypes.size();
  unsigned numHruIdNames = hruIds.size();

  //check if names are consistent
  //for which hrus we want to change the stor type - vector of character
  if(numINTRTypes!=numHruIdNames) {
    Rcpp::Rcout << "The number of INTERCEPTION types does not correspond to the number of HRUs " << numINTRTypes <<"\n";
    Rcpp::stop("\nWrong size number of interception types.\n");
  } else {
    std::vector<std::string> intcpNameStr = Rcpp::as<std::vector<std::string> >(intcptnTypes);
   for(unsigned it=0; it<numINTRTypes;it++ ){
      if ( std::find(allTnterceptionStorTypeNames.begin(), allTnterceptionStorTypeNames.end(), intcpNameStr[it]) == allTnterceptionStorTypeNames.end()) {
        Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of Interception Type Values.\n");
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
     std::map<std::string, interception_STORtype> s_mapStringToINTCPtype_HRUtype = {
       {"Rutter_Gash", interception_STORtype::Rutter_Gash}
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
     std::vector<std::pair<unsigned,interception_STORtype>> intcptnTypesToLoad;


     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToINTCPtype_HRUtype[intcpNameStr[id]]) {
         case interception_STORtype::Rutter_Gash:
           intcptnTypesToLoad.push_back(std::make_pair(indexHru[id], interception_STORtype::Rutter_Gash));
           break;
         }
       }
    dHRUM_ptr.get()->initIntrcptnStypeToAlldHrus(intcptnTypesToLoad);
  }



  return;
}
