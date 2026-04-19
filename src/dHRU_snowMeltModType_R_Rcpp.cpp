#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of snow melt models types to dHRU model for all single HRUs.
//'
//' Setting the snow melt type to dHRUM to all HRUs. Possibe types: \code{DFF}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param snowMeltModelTypes a charater vector of Surface retention type names
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes=rep("DDF",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
 void setSnowMeltModeltypeToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector snowMeltModelTypes, Rcpp::CharacterVector hruIds) {
   unsigned numSnowMeltModelTypes = snowMeltModelTypes.size();
   unsigned numHruIdNames = hruIds.size();

   //check if names are consistent
   //for which hrus we want to change the stor type - vector of character
   if(numSnowMeltModelTypes!=numHruIdNames) {
     Rcpp::Rcout << "The number of Snow Melt Model types does not correspond to the number of HRUs " << numSnowMeltModelTypes <<"\n";
     Rcpp::stop("\nWrong size number of Surface Retention types.\n");
   } else {

     std::vector<std::string> snowMeltModelsnameStr = Rcpp::as<std::vector<std::string> >(snowMeltModelTypes);

     for(unsigned it=0; it<numSnowMeltModelTypes;it++ ){
       if ( std::find(allSnowMdls.begin(), allSnowMdls.end(), snowMeltModelsnameStr[it]) == allSnowMdls.end()) {
         Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
         Rcpp::stop("\n Wrong names of Snow Melt model Type Values.\n");
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

     std::map<std::string, snow_Model> s_mapStringToSMtype_HRUtype = {
       {"DDF", snow_Model::DDF}
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

     std::vector<std::pair<unsigned,snow_Model>> snowMeltModelTypesToLoad;

     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToSMtype_HRUtype[snowMeltModelsnameStr[id]]) {
        case snow_Model::DDF:
         snowMeltModelTypesToLoad.push_back(std::make_pair(indexHru[id], snow_Model::DDF));
         break;
       }
     }
     dHRUM_ptr.get()->initSnowMelMdltypeToAlldHrus(snowMeltModelTypesToLoad);
   }

   return;
 }
