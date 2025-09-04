#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of pond model types to dHRU model for all single HRUs.
//'
//' Setting of params to dHRUM.
//'
//' @param dHRUM_ptr pointer to dHRUM instance
//' @param PondTypes vector of values of parameters
//' @param hruIds a charater vector of parameter names
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
 void setPondToAlldHrus(Rcpp::XPtr<dHRUM> dHRUM_ptr, Rcpp::CharacterVector PondTypes, Rcpp::CharacterVector hruIds) {
   unsigned numPondTypes = PondTypes.size();
   unsigned numHruIdNames = hruIds.size();

   //check if names are consistent
   //for which hrus we want to change the fast response type - vector of character
   if(numPondTypes!=numHruIdNames) {
     Rcpp::Rcout << "The number of ponds does not correspond to the number of HRUs " << numPondTypes <<"\n";
     Rcpp::stop("\nWrong size number of ponds.\n");
   } else {

     std::vector<std::string> PondNameStr = Rcpp::as< std::vector<std::string>>(PondTypes);

     for(unsigned it=0; it<numPondTypes;it++ ){
       if ( std::find(all_pondNames.begin(), all_pondNames.end(), PondNameStr[it]) == all_pondNames.end()) {
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

     std::map<std::string, pond_type> s_mapStringToPondtype_HRUtype = {
       {"noPond", pond_type::noPond},
       {"pondBasic", pond_type::pondBasic},
       {"pondSoilSois", pond_type::pondSoilSois},
       {"pondGWGros", pond_type::pondGWGros}

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

     std::vector<std::pair<unsigned,pond_type>> PondTypesToLoad;

     for(unsigned id=0;id<numHruIdNames;id++) {
       switch(s_mapStringToPondtype_HRUtype[PondNameStr[id]]) {
        case pond_type::noPond:
         PondTypesToLoad.push_back(std::make_pair(indexHru[id], pond_type::noPond));
         break;
       case pond_type::pondBasic:
         PondTypesToLoad.push_back(std::make_pair(indexHru[id], pond_type::pondBasic));
         break;
       case pond_type::pondSoilSois:
         PondTypesToLoad.push_back(std::make_pair(indexHru[id], pond_type::pondSoilSois));
         break;
       case pond_type::pondGWGros:
         PondTypesToLoad.push_back(std::make_pair(indexHru[id], pond_type::pondGWGros));
         break;
       }
     }
     dHRUM_ptr.get()->initPondToAlldHrus(PondTypesToLoad);
   }

   return;
 }
