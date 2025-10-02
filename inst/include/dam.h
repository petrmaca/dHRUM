#ifndef RESERVOIR_H
#define RESERVOIR_H

#include <fstream>
#include <iostream>
#include <iomanip>
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
#include "dta_dam_1d.h"
#include "damSel.h"

class dam {
public:
  dam();
  ~dam();
  dam(const dam& other);
  dam& operator=(const dam& other);

  void set_calender();
  numberSel get_initState(const init_dStype& _Stype);
  void set_ZeroinitStates(const unsigned& numres);
  void init_inputs(numberSel val, unsigned numDTA);
  void set_data(const hdata& dta,const dam_ts& _tsType);


protected:

private:
  numberDta tstRM;//!< The counter for main loop in run model
  dta_dam_1d hyd_dta;//!< Data of all time series variables

  numberSel prev_DamS;//!< The helper variable for updating reservoir storage

  numberSel damMax;//!< Minimum volume of the dam [m3]
  numberSel MRF; //!< Minimum Residual Flow [m3/s]
  numberSel damArea; //!< Area of the dam in [m2] - for groundwater communication, WS evaporation,...
  numberSel damLeng; //!< Length of the dam body[m] - for dam body leakage
  numberSel damBank; //!< Length of the dam bank[m] - for Soil communication

  ETdam_type damEt_TYPE;//!< Type of evaporation
  DamSOISPerc_type damSperc_TYPE;//!< Type of soil percolation
  DamGWPerc_type damGperc_TYPE;//!< Type of groundwater percolation
  DamRouT_type damRout_TYPE;//!< Type of regular outflow



};

#endif // RESERVOIR_H
