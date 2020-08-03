#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <valarray>
#include <vector>
#include <utility> // pair
#include <string> //string

#include "numberSel.h"

class params {
 public:
  params();
  virtual ~params();
  params(const params& other);
  params& operator=(const params& other);

  void s_params(const numberSel& dta,par_HRUtype _parType);//!< The setting of model parameters
  numberSel g_par(const par_HRUtype& _parType);//!< The getting of model parameters

  unsigned g_numFastRes();//!< Get the number of fast reservoirs
  void s_numFastRes(const unsigned& numRes);//!< set number of fast reservoirs


  void s_parLoadToCalib(std::vector<std::pair<numberSel,par_HRUtype>>& parsToLoad);
  void s_params(const std::pair <numberSel, par_HRUtype>& parDta);
  void s_default();

  void p_param();

  unsigned g_numPars();//!< Get the number of parameters

 protected:

 private:
  //!< The PDM parametrization based on the Pareto distribution
  unsigned numPars;//!<  The number of model parameters
  pdata pars;//!< The values of parameters on given HRU
  pdata up_pars;//!< The upper bounds of parameters
  pdata low_pars;//!< The lower bound of parameters
  unsigned numFastRes;//!< The number of fast runoff reservoirs, states are implemented in data_HB_1d class

};

#endif // PARAMS_H
