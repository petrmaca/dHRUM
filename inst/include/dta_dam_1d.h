#ifndef DTA_DAM_1D_H
#define DTA_DAM_1D_H

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

#include "damSel.h"
#include "numberSel.h"

class dta_dam_1d {
public:
  dta_dam_1d();
  dta_dam_1d(unsigned _ndata);
  virtual ~dta_dam_1d();
  dta_dam_1d(const dta_dam_1d& other);
  dta_dam_1d& operator=(const dta_dam_1d& rhs);


  numberSel get_daysInMonth(const unsigned& tstMonth, const unsigned& year);
  bool leap_Check_Year(unsigned TestedYear);
  void s_Julian_day();
  void s_calender();
  // void p_calender();
  void s_initDate(const unsigned& Year, const unsigned& Month,const unsigned& Day,const unsigned& initNumTS);

  void s_initStates(const hdata& initfastRes, const numberSel& init_State,const init_dStype& _Stype);
  numberSel g_initState(const init_dStype& _Stype);

  void s_varVal(const numberSel& dta, const unsigned& tst,const dam_ts& _tsType);
  void s_data(const hdata& dta,const dam_ts& _tsType, bool updateNumTS);
  numberSel g_dta(const unsigned& tst,const dam_ts& _tsType);
  hdata get_HbTsData(const dam_ts& _tsType);
  void setOneTstoZero(const dam_ts& _tsType);

protected:

private:
  unsigned numTS;//!< The number of time intervals

  caldata year;//!< Year
  caldata month;//!< Month
  caldata day;//!< Day
  caldata Jday;//!< Julian day with no translation to CE date in days
  caldata init_year;
  caldata init_month;
  caldata init_day;

  //inflows
  hdata Prec;//!< Precipitation on water surface
  hdata InfL;//!< Inflow from water channels
  hdata InlT;//!< Water inlet (transfer between two dams)[m3/day]

  //outflows
  hdata EtdM;//!< Evaporation from water surface [mm/day]
  hdata OufL;//!< Total outflow
  hdata OflW;//!< overflow [m3/day]
  hdata OulT;//!< Water outlet (transfer between dam and other dam or industry) [m3/day]

  //reversible flows
  hdata DaiS;//!< Soil percolation [m/s]
  hdata DaiG;//!< Groundwater percolation [m/s]

  //Storage
  hdata DamS;//!< Dam storage[m3]
  numberSel init_DamS;//!< Initial value of dam storage

};

#endif // DTA_DAM_1D_H
