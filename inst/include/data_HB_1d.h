#ifndef DATA_HB_1D_H
#define DATA_HB_1D_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <valarray>
#include <cmath>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include "numberSel.h"

class data_HB_1d {
 public:
  data_HB_1d();
  data_HB_1d(unsigned _ndata);
  virtual ~data_HB_1d();
  data_HB_1d(const data_HB_1d& other);
  data_HB_1d& operator=(const data_HB_1d& rhs);

  void s_data(const hdata& dta,const ts_type& _tsType, bool updateNumTS);//!< Setter  of private vectors of RR data
  void s_varVal(const numberSel& dta, const unsigned& tst,const ts_type& _tsType);

  void s_initStates(const hdata& initfastRes, const numberSel& init_State,const init_Stype& _Stype);
  void s_init_states_fastRes(const hdata& initfastRes);

  numberSel g_dta(const unsigned& tst,const ts_type& _tsType);
  numberSel g_initState(const init_Stype& _Stype);
  hdata g_initFastresStates();
  numberSel g_oneFastResOut(const unsigned& itFasRes);
  void s_oneFastreResUout(const numberSel& resOut,const unsigned& itFasRes);
  numberSel g_oneFastResState(const unsigned& itFasRes);
  void s_oneFastreResUstate(const numberSel& resState,const unsigned& itFasRes);

  unsigned g_numdta();

  //!< The PET calculations
  void s_Pet_Pars(const numberSel& newLatitude, const pet_Type& newPeType);
  void calc_Pet();
  void OudinPET();
  void HamonPET();
  bool leap_Check_Year(unsigned TestedYear);
  void s_Julian_day();

  void s_calender();
  void p_calender();
  void s_initDate(const unsigned& Year, const unsigned& Month,const unsigned& Day,const unsigned& initNumTS);

  numberDta g_calDta(const cal_Type& calDate, const unsigned& ts);

  hdata get_HbTsData(const ts_type& _tsType);
  caldata getCalData(const cal_Type& calDate);
  void setAllToZeros(const bool& CalDta,const bool& TsDta);
  void setOneTstoZero(const ts_type& _tsType);
  void setOneCalDateToZero(const cal_Type& calDate);

  void printDataToFile(const std::string& Filet);

  void loadCalData(const caldata& yyear, const caldata& mmonth, const caldata& dday);

 protected:

 private:
  unsigned numTS;//!< The number of time intervals

  numberSel Latitude;//!< The latitude og basin centroid in degrees
  pet_Type PETtype;//!< The tag on selected PET method

  caldata year;//!< Year
  caldata month;//!< Month
  caldata day;//!< Day
  caldata Jday;//!< Julian day with no translation to CE date in days
  caldata init_year;
  caldata init_month;
  caldata init_day;

  hdata Prec;//!< Precipitation
  hdata Snow;//!< Snow depth
  hdata AEt;//!< Actual Evapotranspiration
  hdata PEt;//!< Potential Evapotranspiration

  hdata Temp;//!< Temperature

  hdata TroF;//!< Throughfall
  hdata SteF;//!< Stemflow
  hdata CanF;//!< Canopy drainage

  hdata CanS;//!< Canopy storage
  hdata SteS;//!< Stem storage

  hdata EvaC;//!< Canopy Evaporation
  hdata EvaS;//!< Stem Evaporation
  hdata EvbS;//!< Bare soil Evapotranspiration
  hdata IntS;//!< Interception storage

  hdata SoiS;//!< Soil storage
  hdata GroS;//!< Groundwater storage
  hdata GroS1;//!< Groundwater storage 1
  hdata GroS2;//!< Groundwater storage 2
  hdata SurS;//!< Surface retention

  hdata TotR;//!< Total runoff
  hdata Basf;//!< Baseflow
  hdata DirR;//!< Direct Runoff

  hdata Melt;//!< SnowMelt
  hdata Perc;//!< Percolation
  hdata Pref;//!< Effective Precipitation

  numberSel init_SoiS;//!< Initial value of soil storage
  numberSel init_GroS;//!< Initial value of groundwater storage
  numberSel init_GroS1;//!< Initial value of groundwater storage 1
  numberSel init_GroS2;//!< Initial value of groundwater storage 2
  numberSel init_CanS;//!< Initial value of Canopy Interception storage
  numberSel init_SteS;//!< Initial value of Stem Interception storage
  numberSel init_SnoS;//!< Initial variable of Snow storage
  numberSel init_SurS;//!< Initial value of Surface retention storage

  unsigned numfastRes;//!<  the number of fast response reservoirs
  hdata StateFastRes;//!< The state variables of fast response reservoirs
  hdata OutFastRes;//!< The outputs of serie for fast response reservoirs
};

#endif // DATA_HB_1D_H
