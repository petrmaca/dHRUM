#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of surface retention models types to dHRU model for all single HRUs.
//'
//' Setting the fast response type to dHRUM to all HRUs. Possibe types: \code{SerialCascadeLinRes}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param fastResponseTypes a character vector of fast responses
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
 void setFastResponsesToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector fastResponseTypes, Rcpp::CharacterVector hruIds) {
   unsigned numFASTRESPONSESTYPES = fastResponseTypes.size();
   unsigned numHruIdNames = hruIds.size();

   //check if names are consistent
   //for which hrus we want to change the fast response type - vector of character
   if(numFASTRESPONSESTYPES!=numHruIdNames) {
     Rcpp::Rcout << "The number of Surface Retention types does not correspond to the number of HRUs " << numFASTRESPONSESTYPES <<"\n";
     Rcpp::stop("\nWrong size number of Surface Retention types.\n");
   } else {

     std::vector<std::string> fastResponsesNameStr = Rcpp::as< std::vector<std::string>>(fastResponseTypes);

     for(unsigned it=0; it<numFASTRESPONSESTYPES;it++ ){
       if ( std::find(all_FastResponsesNames.begin(), all_FastResponsesNames.end(), fastResponsesNameStr[it]) == all_FastResponsesNames.end()) {
         Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
         Rcpp::stop("\n Wrong names of Surface Retention Type Values.\n");
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

     std::map<std::string, fast_Response> s_mapStringToINTCPtype_HRUtype = {
       {"SerialCascadeLinRes", fast_Response::SerialCascadeLinRes},
       {"SerialLinResGWPerc", fast_Response::SerialLinResGWPerc}
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

     std::vector<std::pair<unsigned,fast_Response>> fastResponseTypesToLoad;

     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToINTCPtype_HRUtype[fastResponsesNameStr[id]]) {
        case fast_Response::SerialCascadeLinRes:
         fastResponseTypesToLoad.push_back(std::make_pair(indexHru[id], fast_Response::SerialCascadeLinRes));
         break;
       case fast_Response::SerialLinResGWPerc:
         fastResponseTypesToLoad.push_back(std::make_pair(indexHru[id], fast_Response::SerialLinResGWPerc));
         break;
       }
     }
     dHRUM_ptr.get()->initFastResponsesToAlldHrus(fastResponseTypesToLoad);
   }

   return;
 }
