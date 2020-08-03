#ifndef DHRUM_H
#define DHRUM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <map>
#include <omp.h>
// #include <cmath>

#include "numberSel.h"
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
  void setInputsToAllHrus(std::string namesFilePath);
  void loadPTDatToAllHrus(hdata Prec, hdata Temp, const numberSel& val,const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay);
  void setInputsToOneHru(std::string namesFilePath, unsigned Id);
  void setParamsToAllHrus(std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad);
  void setParamsToOneHru(std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad, unsigned Id);
  void calcPetToAllHrus(numberSel Latit, pet_Type myPetType = pet_Type::OUDIN);
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

protected:

private:
  std::vector<single_HMunit> dHruVec;//!< vector of single HRU objects
  numberDta dimHM;//!< Number of single HRU objects
  numberSel basinArea;//!< The total basin area
  hdata Areas;//!< Arrays of Areas of all single HRUs
  numberDta numTs;//!< Number of time intervals

  data_HB_1d basinDta;//!< The HB data on all basin

};

#endif // DHRUM_H
