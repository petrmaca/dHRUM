#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of interception models types to dHRU model for all single HRUs.
//'
//' Setting the interception type to dHRUM to all HRUs. Possible types: \code{Rutter_Gash,van_Dijk,Eliades}
//' Setting the model for calculating the Smac based on LAI models: \code{Pitman,VonHoyningenHuene}
//'
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param intcptnTypes a charater vector of Interception type names
//' @param hruIds ids on Hrus
//' @param InstStLai the TRUE/FALSE vector allowing the use of the LAI trnasformed max interception, canopy, and stem storage
//' @param smaxlaiTypes the name of the model used for calculating the Smax based on lai
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
//' smaxlaiTypes = rep("VonHoyningenHuene",times= length(Areas))
//' InstStLai = rep(TRUE,times= length(Areas))
//' setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus,InstStLai,smaxlaiTypes )
//'
//' setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus)
// [[Rcpp::export]]
void setInterceptiontypeToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector intcptnTypes, Rcpp::CharacterVector hruIds, Rcpp::LogicalVector InstStLai, Rcpp::CharacterVector smaxlaiTypes) {
  unsigned numINTRTypes = intcptnTypes.size();
  unsigned numHruIdNames = hruIds.size();
  unsigned numIntStLai = InstStLai.size();

  //check if names are consistent
  //for which hrus we want to change the stor type - vector of character
  if((numINTRTypes!=numHruIdNames)||(numINTRTypes!=numIntStLai)||(numHruIdNames!=numIntStLai) ){
    Rcpp::Rcout << "The number of INTERCEPTION types does not correspond to the number of HRUs or number of Lai TRUE/FALSE" << numINTRTypes <<"\n";
    Rcpp::stop("\nWrong size number of interception types or LAI settings.\n");
  } else {
    std::vector<std::string> intcpNameStr = Rcpp::as<std::vector<std::string> >(intcptnTypes);

    for(unsigned it=0; it<numINTRTypes;it++ ){
      if ( std::find(allTnterceptionStorTypeNames.begin(), allTnterceptionStorTypeNames.end(), intcpNameStr[it]) == allTnterceptionStorTypeNames.end()) {
        Rcpp::Rcout << "\nSomething wrong on item " << (it+1) << "\n";
        Rcpp::stop("\n Wrong names of Interception Type Values.\n");
      }
    }

    for(unsigned it = 0; it < numIntStLai; ++it) {

      if (Rcpp::LogicalVector::is_na(InstStLai[it])) {
        Rcpp::Rcout << "TRUE/FALSE on LAI for HRU " << it << " is NA\n";
        Rcpp::stop("\nWrong inputs of interception types or LAI settings.\n");
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
       {"Rutter_Gash", interception_STORtype::Rutter_Gash},
       {"van_Dijk", interception_STORtype::van_Dijk},
       {"Eliades", interception_STORtype::Eliades}
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
       case interception_STORtype::van_Dijk:
         intcptnTypesToLoad.push_back(std::make_pair(indexHru[id], interception_STORtype::van_Dijk));
         break;
       case interception_STORtype::Eliades:
         intcptnTypesToLoad.push_back(std::make_pair(indexHru[id], interception_STORtype::Eliades));
         break;
         }
       }

     std::vector<bool> vecINtStLai;
     vecINtStLai.resize(numIntStLai);

     for(unsigned it = 0; it < numIntStLai; ++it) {
       if (Rcpp::LogicalVector::is_na(InstStLai[it])) {
         vecINtStLai[it] = false; // Custom rule: Treat NA as false
       } else {
         vecINtStLai[it] = InstStLai[it];
       }
     }

     std::map<std::string, lai_SmaxModel> s_mapStringToSmaxLai_HRUtype = {
       {"Pitman", lai_SmaxModel::Pitman},
       {"VonHoyningenHuene", lai_SmaxModel::VonHoyningenHuene}
     };


     std::vector<std::pair<unsigned,lai_SmaxModel>> SmaxLaiTypesToLoad;

     std::vector<std::string> smaxlaiNameStr = Rcpp::as<std::vector<std::string> >(smaxlaiTypes);

     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToSmaxLai_HRUtype[smaxlaiNameStr[id]]) {
       case lai_SmaxModel::Pitman:
         SmaxLaiTypesToLoad.push_back(std::make_pair(indexHru[id], lai_SmaxModel::Pitman));
         break;
       case lai_SmaxModel::VonHoyningenHuene:
         SmaxLaiTypesToLoad.push_back(std::make_pair(indexHru[id],  lai_SmaxModel::VonHoyningenHuene));
         break;
       }
     }


    dHRUM_ptr.get()->initIntrcptnStypeToAlldHrus(intcptnTypesToLoad, vecINtStLai,SmaxLaiTypesToLoad);
  }



  return;
}
