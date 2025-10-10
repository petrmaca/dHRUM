#include <Rcpp.h>
#include "dHRUM.h"
//' Sets the types of pond model types to dHRU model for all single HRUs.
//'
//' Setting of pond to dHRUM.
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
//' setPondToAlldHrus(dHRUM_ptr = dhrus,PondTypes=rep("Pond",times= length(Areas)),hruIds=IdsHrus)
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
       {"Pond", pond_type::Pond}

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
       case pond_type::Pond:
         PondTypesToLoad.push_back(std::make_pair(indexHru[id], pond_type::Pond));
         break;
       }
     }
     dHRUM_ptr.get()->initPondToAlldHrus(PondTypesToLoad);
   }

   return;
 }




//' Sets the types of pond to dHRU model for one single HRU.
//'
//' Setting of pond to one particular single HRU for dHRU.
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
//' DateVec = as.Date(c("1990/01/30","1990/01/31","1990/02/01")))
//' ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
//' SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
//' pondDF1 = data.frame( pondArea = 40500, PonsMax= 45000, MRF= 0.039)
//' pondDF2 = data.frame( ET = "ETpond1", in_SOISperc= "noPondSOISPerc", in_GWperc= "noPondGWPerc", out_SOISperc= "noPondSOISPerc", out_GWperc= "noPondGWPerc",regular_out="PondRouT3" )
//' HruPondID = 1
//' setPondToOnedHru(dHRUM_ptr = dhrus,HruPondID,names(pondDF1),as.numeric(pondDF1),as.character(pondDF2),names(pondDF2))
// [[Rcpp::export]]
void setPondToOnedHru(Rcpp::XPtr<dHRUM> dHRUM_ptr,unsigned singleHruId, Rcpp::CharacterVector ValNames,Rcpp::NumericVector ValVals,Rcpp::CharacterVector TypeNames,Rcpp::CharacterVector TypeVals) {


   unsigned numValNames = ValNames.size();
   unsigned numValVals = ValVals.size();
   unsigned numTypeNames = TypeNames.size();
   unsigned numTypeVals = TypeVals.size();
   unsigned numPars1Hru;
   numPars1Hru = dHRUM_ptr.get()->get_singleHRUnumPars(singleHruId);

    if((numValVals+numTypeVals)!=9){
     Rcpp::Rcout << "The number of set pond inputs is: " << (numValVals+numTypeVals) <<"\n";
     Rcpp::stop("\n but  required number of pond inputs is: 9.\n");
   } else  {

     std::vector<numberSel> valvals = Rcpp::as<std::vector<numberSel>>(ValVals);
     std::vector<std::string>  valnames = Rcpp::as<std::vector<std::string> >(ValNames);

     std::vector<std::pair<std::string,numberSel>> PondDefs;//pond definitions
     for(unsigned it=0; it<numValNames;it++ ){
       PondDefs.push_back(std::make_pair(valnames[it],valvals[it]));
     }

     std::vector<std::string> typevals = Rcpp::as<std::vector<std::string>>(TypeVals);
     std::vector<std::string> typenames = Rcpp::as<std::vector<std::string> >(TypeNames);

     std::vector<std::pair<std::string,std::string>> PondBeh;//pond behavior
     for(unsigned it=0; it<numTypeNames;it++ ){
       PondBeh.push_back(std::make_pair(typevals[it],typenames[it]));
     }


     unsigned dHRUdim = 0;
     dHRUdim = dHRUM_ptr.get()->getdHRUdim();
     if( singleHruId > dHRUdim){
       Rcpp::stop("The wrong ID or number of single HRU.\n");
     }

     dHRUM_ptr.get()->initPondToOneHRU(singleHruId,PondDefs,PondBeh);
   }
   return  ;
 }
