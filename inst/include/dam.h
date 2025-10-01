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
#include "damSel.h"

class dam {
public:
  dam();
  ~dam();
  dam(const dam& other);
  dam& operator=(const dam& other);

  void set_calender();

protected:

private:
  dta_dam_1d hyd_dta;//!< The data of all time series of hydrological variables

  numberSel prev_damStor;//!< The helper variable for updating reservoir storage
  numberSel damArea;//!< The area of reservoir unit in m2



};

#endif // RESERVOIR_H
