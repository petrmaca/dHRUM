#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of surface retention models types to dHRU model for all single HRUs.
//'
//' Setting the surface retention type to dHRUM to all HRUs. Possibe types: \code{SurfaceAll}, \code{SurfacePRTL}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param surfaceStorTypes a charater vector of Surface retention type names
//' @param hruIds ids on Hrus
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
 void setSurfaceStortypeToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector surfaceStorTypes, Rcpp::CharacterVector hruIds) {
   unsigned numSURFSTOTypes = surfaceStorTypes.size();
   unsigned numHruIdNames = hruIds.size();

   //check if names are consistent
   //for which hrus we want to change the stor type - vector of character
   if(numSURFSTOTypes!=numHruIdNames) {
     Rcpp::Rcout << "The number of Surface Retention types does not correspond to the number of HRUs " << numSURFSTOTypes <<"\n";
     Rcpp::stop("\nWrong size number of Surface Retention types.\n");
   } else {

     std::vector<std::string> surfaceStorNameStr = Rcpp::as<std::vector<std::string> >(surfaceStorTypes);

     for(unsigned it=0; it<numSURFSTOTypes;it++ ){
       if ( std::find(allSurfacStorsNames.begin(), allSurfacStorsNames.end(), surfaceStorNameStr[it]) == allSurfacStorsNames.end()) {
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

     std::map<std::string, surface_STORtype> s_mapStringToINTCPtype_HRUtype = {
       {"SurfaceAll", surface_STORtype::SurfaceAll},
       {"SurfacePRTL", surface_STORtype::SurfacePRTL},
       {"Wetland", surface_STORtype::Wetland},
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

     std::vector<std::pair<unsigned,surface_STORtype>> surfaceStorTypesToLoad;

     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToINTCPtype_HRUtype[surfaceStorNameStr[id]]) {
        case surface_STORtype::SurfaceAll:
         surfaceStorTypesToLoad.push_back(std::make_pair(indexHru[id], surface_STORtype::SurfaceAll));
         break;
        case surface_STORtype::SurfacePRTL:
         surfaceStorTypesToLoad.push_back(std::make_pair(indexHru[id], surface_STORtype::SurfacePRTL));
         break;
       case surface_STORtype::Wetland:
         surfaceStorTypesToLoad.push_back(std::make_pair(indexHru[id], surface_STORtype::Wetland));
         break;
       }
     }
     dHRUM_ptr.get()->initSurfaceStypeToAlldHrus(surfaceStorTypesToLoad);
   }

   return;
 }
