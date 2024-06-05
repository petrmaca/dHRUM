#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the similar values of params to dHRU model for all single HRUs.
//'
//' Setting the groudnwater type of dHRUM equal to all HRUs. Possibe types: \code{LIN_RES} \code{LINL_RES} \code{LINBY_RES}
//' \code{POW_RES}, \code{EXP_RES} \code{LIN_2SE} \code{LIN_2PA} \code{FLEX_RES} \code{EXP_LOG}
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param gwTypes a charater vector of GW type names
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
void setGWtypeToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector gwTypes, Rcpp::CharacterVector hruIds) {
  unsigned numGWTypes = gwTypes.size();
  unsigned numHruIdNames = hruIds.size();

  //check if names are consistent
  //for which hrus we want to change the stor type - vector of character
  if(numGWTypes!=numHruIdNames) {
    Rcpp::Rcout << "The number of GW types does not correspond to the number of HRUs " << numGWTypes <<"\n";
    Rcpp::stop("\nWrong size number of gw types.\n");
  } else {
    std::vector<std::string> gwNameStr = Rcpp::as<std::vector<std::string> >(gwTypes);
   for(unsigned it=0; it<numGWTypes;it++ ){
      if ( std::find(allGWStorTypeNames.begin(), allGWStorTypeNames.end(), gwNameStr[it]) == allGWStorTypeNames.end()) {
        Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of GW Type Values.\n");
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
     std::map<std::string, gs_STORtype> s_mapStringToGWtype_HRUtype = {
       {"LIN_RES", gs_STORtype::LIN_RES},
       {"LINL_RES", gs_STORtype::LINL_RES },
       {"LINBY_RES",  gs_STORtype::LINBY_RES},
       {"POW_RES",  gs_STORtype::POW_RES},
       {"EXP_RES",  gs_STORtype::EXP_RES},
       {"LIN_2SE",  gs_STORtype::LIN_2SE},
       {"LIN_2PA",  gs_STORtype::LIN_2PA},
       {"FLEX_RES", gs_STORtype::FLEX_RES},
       {"EXP_LOG", gs_STORtype::EXP_LOG}
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
     std::vector<std::pair<unsigned,gs_STORtype>> gwTypesToLoad;


     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToGWtype_HRUtype[gwNameStr[id]]) {
         case gs_STORtype::LIN_RES:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::LIN_RES));
           break;
         case gs_STORtype::LINL_RES:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::LINL_RES));
           break;
         case gs_STORtype::LINBY_RES:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::LINBY_RES));
           break;
         case gs_STORtype::POW_RES:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::POW_RES));
           break;
         case gs_STORtype::EXP_RES:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::EXP_RES));
           break;
         case gs_STORtype::LIN_2SE:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::LIN_2SE));
           break;
         case gs_STORtype::LIN_2PA:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::LIN_2PA));
           break;
         case gs_STORtype::FLEX_RES:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::FLEX_RES));
           break;
         case gs_STORtype::EXP_LOG:
           gwTypesToLoad.push_back(std::make_pair(indexHru[id], gs_STORtype::EXP_LOG));
           break;
         }
       }
    dHRUM_ptr.get()->initGWtypeToAlldHrus(gwTypesToLoad);
  }



  return;
}
