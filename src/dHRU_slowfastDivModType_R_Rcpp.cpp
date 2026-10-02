#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of snow melt models types to dHRU model for all single HRUs.
//'
//' Setting the snow melt type to dHRUM to all HRUs. Possibe types: \code{DFF}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param adivsModelTypes a charater vector of Surface retention type names
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setFastSlowDivModeltypeToAlldHrus(dHRUM_ptr = dhrus,adivsModelTypes=rep("CnstAdiv",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
 void setFastSlowDivModeltypeToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector adivsModelTypes, Rcpp::CharacterVector hruIds) {
   unsigned numadivsModelTypes = adivsModelTypes.size();
   // std::cout << numadivsModelTypes << std::endl;
   unsigned numHruIdNames = hruIds.size();

   //check if names are consistent
   //for which hrus we want to change the stor type - vector of character
   if(numadivsModelTypes!=numHruIdNames) {
     Rcpp::Rcout << "The number of Adiv models of divider of slow fast responce Model types does not correspond to the number of HRUs " << numadivsModelTypes <<"\n";
     Rcpp::stop("\nWrong size number of slow fast responce divider model types.\n");
   } else {

     std::vector<std::string> adivModelsnameStr = Rcpp::as<std::vector<std::string> >(adivsModelTypes);

     for(unsigned it=0; it<numadivsModelTypes;it++ ){
       if ( std::find(allAdivMdls.begin(), allAdivMdls.end(), adivModelsnameStr[it]) == allAdivMdls.end()) {
         Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
         Rcpp::stop("\n Wrong names of adiv slow/fast model Type Values.\n");
       }
     }

     const std::vector<std::string> ids = dHRUM_ptr.get()->getHRUIds();
     std::vector<std::string> hruIdName = Rcpp::as<std::vector<std::string> >(hruIds);

     //std::vector<std::string> ids = Rcpp::as<std::vector<std::string> >(hruIds);
     for(unsigned it=0; it<numHruIdNames;it++ ){
       if ( std::find(ids.begin(), ids.end(), hruIdName[it]) == ids.end()) {
         Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
         Rcpp::stop("\n Wrong names of Hru Id Values when looking at setter of adiv types.\n");
       }
     }

     std::map<std::string, adiv_Model> s_mapStringToAdivtype_HRUtype = {
       {"CnstAdiv", adiv_Model::CnstAdiv},
       {"adivSoilSat", adiv_Model::adivSoilSat}
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

     std::vector<std::pair<unsigned,adiv_Model>> adivModelTypesToLoad;

     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToAdivtype_HRUtype[adivModelsnameStr[id]]) {
        case adiv_Model::CnstAdiv:
          adivModelTypesToLoad.push_back(std::make_pair(indexHru[id], adiv_Model::CnstAdiv));
         break;
       case adiv_Model::adivSoilSat:
         adivModelTypesToLoad.push_back(std::make_pair(indexHru[id], adiv_Model::adivSoilSat));
         break;
       }
     }
     dHRUM_ptr.get()->initAdivMdltypeToAlldHrus(adivModelTypesToLoad);
   }

   return;
 }
