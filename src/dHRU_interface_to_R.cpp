#include <Rcpp.h>
#include "dHRUM.h"

//' Initialization of dHRU pointer to a dHRUM
//'
//' Creates pointer instance of dHRUM for the catchment.
//' initializes a dimension of dHRUM controled by the number of single Hru, areas of all single Hrus,
//' and ID's of all single Hrus.
//'
//' @param dimdHru a single \code{numberDta} number of single Hrus.
//' @param vecAreas a \code{numeric vector} of size \code{dimHru} of Areas for all single HRUs on dHRU.
//' @param hrusIDs a \code{character vector} of size \code{dimHru} of Id's for all single HRUs on dHRU.
//' @param nthreads the number of r dHRUM instance
//' @return dHRUM_ptr pointer to dHru instance.
//' @export
//' @examples
//' nHrus <- 200
//' Areas <- runif(nHrus,min = 1,max  = 10)
//' IdsHrus <- paste0("ID",seq(1:length(Areas)))
//' dhrus <- initdHruModel(nHrus,Areas,IdsHrus,1)
// [[Rcpp::export]]
Rcpp::XPtr<dHRUM> initdHruModel(numberDta dimdHru, Rcpp::NumericVector vecAreas, Rcpp::StringVector hrusIDs,unsigned nthreads) {
  // implicit form
  // 1) creates dHRU instance on the heap (allocates memory and call constructor with no arguments)
  // 2) creates dhruPTR variable on the stack initialized with pointer to dHRU's instance
  dHRUM* dHRUM_ptr = new dHRUM();

  unsigned vecArSize = vecAreas.size();
  unsigned vecIdNamesSize = hrusIDs.size();
  if(!((dimdHru == vecArSize) && (dimdHru == vecIdNamesSize) && (vecIdNamesSize == vecArSize))) {
    Rcpp::stop("\nThe dim of dHRU single HRU units does not correspond to the length of Areas or Ids.\n");
  } else {
    dHRUM_ptr->set_num_treads(nthreads);
    single_HMunit sHRU_to_VEC;
    dHRUM_ptr->initHrusVec(dimdHru, sHRU_to_VEC);
    std::vector<std::string> vecIDs;
    vecIDs.resize(dimdHru);
    for(unsigned id=0;id<dimdHru;id++){
      vecIDs[id] =  hrusIDs[id];
    }
    hdata vecAreasHD(1,dimdHru);
    for(unsigned aa=0;aa<dimdHru;aa++){
      vecAreasHD[aa] = (numberSel) vecAreas[aa];
      // Rcpp::Rcout << vecAreasHD[aa] << "\n";
    }
    dHRUM_ptr->initHrusID(vecIDs);
    // for(unsigned aa=0;aa<dimdHru;aa++){
    //   Rcpp::Rcout << "Id of  single Hru on positon " << aa << " is  " << dHRU_ptr->getSingleHruId(aa) << "\n";
    // }
    dHRUM_ptr->setAreasToHrus(vecAreasHD);
    // Rcpp::Rcout << "The number of threads is set to " << dHRUM_ptr->get_num_treads() << "\n";
  }
  return Rcpp::XPtr<dHRUM>(dHRUM_ptr);
}
