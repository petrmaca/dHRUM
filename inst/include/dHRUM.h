#ifndef DHRUM_H
#define DHRUM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <string>
#include <map>
#include <omp.h>
// #include <cmath>

#include "numberSel.h"
//#include "pondSel.h"
#include "single_HMunit.h"


class dHRUM {
 public:
  dHRUM();
  ~dHRUM();
  dHRUM(const dHRUM& other);
  dHRUM& operator=(const dHRUM& other);

  void initHrusVec(const unsigned& numHMunits, const single_HMunit& shru);
  void initHrusID(const std::vector<std::string>& SingleHrusIds);
  std::string getSingleHruId(unsigned hruId);
  std::vector<std::string> getHRUIds();
  void setInputsToAllHrus(std::string namesFilePath);
  void loadPTDatToAllHrus(hdata Prec, hdata Temp, const numberSel& val,const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay);
  void setInputsToOneHru(std::string namesFilePath, unsigned Id);
  void loadPTInputsToOneHru(hdata Prec, hdata Temp, const numberSel& val,const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay, unsigned HruIt);
  void setParamsToAlldHrus(std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad);
  void setParamsToOneHru(std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad, unsigned Id);
  std::vector<std::string> getRequiredParamsForHru(unsigned Id);

  void calcPetToAllHrus(numberSel Latit, pet_Type myPetType = pet_Type::OUDIN);
  void calcPetToOneHru(numberSel Latit, pet_Type myPetType, unsigned HruId);
  void calcPetToAllHrusDist(hdata LatitVec, std::vector<std::string> petType);

  void calcHbToAllHrus();
  void setAreasToHrus(hdata vec_Areas);
  void setBasinArea();
  void gatherTsFromHrus();
  numberDta get_numTS();
  void set_numTS();
  void initdHRUbasinDTA();

  void printAllDta(const std::string& Filet);
  hdata get_HbDta(const ts_type& _tsType);
  caldata get_CalDta(const cal_Type& _calType);
  numberDta getdHRUdim();
  void load_CalDataToAllHrus(const caldata& yyear, const caldata& mmonth, const caldata& dday);
  void load_PrecTempToAllHrus(const hdata& Prec, const hdata& Temp);
  unsigned get_singleHRUnumPars(unsigned Id);
  void set_numPars();
  void print_Pars();

  void initGWtypeToAlldHrus(std::vector<std::pair<unsigned, gs_STORtype>>& gs_STORtypes);
  void initSoilStypeToAlldHrus(std::vector<std::pair<unsigned,soil_STORtype>>& soil_STORtypes);
  void initIntrcptnStypeToAlldHrus(std::vector<std::pair<unsigned,interception_STORtype>>& interception_STORtype);
  void initSurfaceStypeToAlldHrus(std::vector<std::pair<unsigned,surface_STORtype>>& surface_STORtype);
  void initFastResponsesToAlldHrus(std::vector<std::pair<unsigned,fast_Response>>& fast_RESPONSESTypes);
  void initPondToAlldHrus(std::vector<std::pair<unsigned,pond_type>>& pondTypes);
  void initPondToOneHRU(unsigned hruId,std::vector<std::pair<std::string,numberSel>>& PondDefs,std::vector<std::pair<std::string,std::string>>& PondBeh);

  numberSel getTsDta(const ts_type& _tsType, const unsigned& HruIndex, const unsigned& tst);

  void set_numFastReservoirsToHrus();
  void set_numFastReservoirsVEC(caldata numFR);
  std::vector<gs_STORtype> get_STORtypes();

  unsigned get_num_treads();
  void set_num_treads(const unsigned&  numTHR);

  std::vector<double> get_param_vec(unsigned hruId);
  std::vector<double> get_upparam_vec(unsigned hruId);
  std::vector<double> get_lowparam_vec(unsigned hruId);
  std::vector<std::string> get_param_names(unsigned hruId);
  void Current_Params(unsigned hruId);


protected:

private:
  std::vector<single_HMunit> dHruVec;//!< vector of single HRU objects
  std::vector<std::string> dHruVecId;//!< vector of string HRU ids
  numberDta dimHM;//!< Number of single HRU objects
  numberSel basinArea;//!< The total basin area
  hdata Areas;//!< Arrays of Areas of all single HRUs
  numberDta numTs;//!< Number of time intervals
  caldata numParsAllHRus;//< Vector with Numbers of params

  unsigned num_threads;

  data_HB_1d basinDta;//!< The HB data on all basin
  std::vector<gs_STORtype> gs_STORtypes;//!< The vector on groundwater storage types
  std::vector<soil_STORtype> sw_STORtypes;//!< The vector on soil water storage types
  std::vector<interception_STORtype> interception_STORtypes;//!< The vector on interception storage types
  std::vector<surface_STORtype> surf_STORtypes;//!< The vector on surface retentions type in HRus
  std::vector<fast_Response> fast_RESPONSESTypes;//!< The vector on fast responses types in HRus
  std::vector<pond_type> pondTypes;//!< The vector on fast responses types in HRus

  caldata NumFastRes;//!< The number of fastre reservoirs in serie
};

#endif // DHRUM_H
