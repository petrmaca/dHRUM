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
#include "data_dam_1d.h"

class dam {
public:
  dam();
  ~dam();
  dam(const dam& other);
  dam& operator=(const dam& other);

  numberSel get_dta(const unsigned& tst, const ts_type& _tsType);//!< Getter data on HB of single pdm unit
  void set_data_prec_temp(const hdata& _prec_dta,const hdata& _temp_dta);//!< Setter data on HB of single pdm unit
  void set_data(const hdata& dta,const ts_type&_tsType);//!< Setter data on HB of single pdm unit
  void set_varValue(const numberSel& dta,const  unsigned& tst,const ts_type& _tsType);

  void set_calender();



protected:

private:
  data_dam_1d hyd_dta;//!< The data of all time series of hydrological variables

  numberSel prev_Stor;//!< The helper variable for updating reservoir storage


  numberSel Area;//!< The area of reservoir unit in m2





};

#endif // RESERVOIR_H
