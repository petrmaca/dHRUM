#ifndef DATA_DAM_1D_H
#define DATA_DAM_1D_H

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

class data_dam_1d {
public:
  data_dam_1d();
  data_dam_1d(unsigned _ndata);
  virtual ~data_dam_1d();
  data_dam_1d(const data_dam_1d& other);
  data_dam_1d& operator=(const data_dam_1d& rhs);

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
  void s_numFastRes(const numberDta& nFastRes);
  numberDta g_numFastRes();

  unsigned g_numdta();





  numberSel get_daysInMonth(const unsigned& tstMonth, const unsigned& year);


  bool leap_Check_Year(unsigned TestedYear);
  void s_Julian_day();


  void s_calender();
  // void p_calender();
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
  hdata SoiS;//!< Soil storage
  hdata GroS;//!< Groundwater storage
  hdata TotR;//!< Total runoff

  hdata Etsw;//!< Evaporation from surface water retention
  hdata DamS;//!< Dam storage

};

#endif // DATA_DAM_1D_H
