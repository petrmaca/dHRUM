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


  numberSel get_daysInMonth(const unsigned& tstMonth, const unsigned& year);
  bool leap_Check_Year(unsigned TestedYear);
  void s_Julian_day();
  void s_calender();
  // void p_calender();
  void s_initDate(const unsigned& Year, const unsigned& Month,const unsigned& Day,const unsigned& initNumTS);


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

  //flows
  hdata InfL;//!< Inflow from water channels
  hdata Prec;//!< Precipitation on water surface
  hdata Sois;//!< Soil percolation [m/s]
  hdata Gros;//!< Groundwater percolation [m/s]
  hdata OufL;//!< Total outflow
  hdata EtdM;//!< Evaporation from water surface [mm/day]
  hdata OflW;//!< overflow [m3/day]

  //extra flows
  hdata OulT;//!< Water outlet (transfer between dam and other dam or industry) [m3/day]
  hdata InlT;//!< Water inlet (transfer between two dams)[m3/day]

  //Storage
  hdata DamS;//!< Dam storage[m3]

  //constants
  numberSel damMax;//!< Minimum volume of the dam [m3]
  numberSel MRF; //!< Minimum Residual Flow [m3/s]
  numberSel damArea; //!< Area of the dam in [m2]
  numberSel damLeak; //!< Area of the dam in [m3/s  per 1 m of dam length]

  numberSel init_DamS;//!< Initial value of dam storage

};

#endif // DATA_DAM_1D_H
