#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <valarray>
#include <vector>
#include <utility> // pair
#include <string> //string
#include <list>

#include "numberSel.h"
#include "parStructSels.h"


class params {
 public:
  params();
  virtual ~params();
  params(const params& other);
  params& operator=(const params& other);

  void s_params(const numberSel& dta,par_HRUtype _parType);//!< The setting of model parameters
  numberSel g_par(const par_HRUtype& _parType);//!< The getting of model parameters
  numberSel g_par_up(const par_HRUtype& _parType);//!< Getting upper bounds on model parameters
  numberSel g_par_low(const par_HRUtype& _parType);//!< Getting lower bounds on model parameters

  unsigned g_numFastRes();//!< Get the number of fast reservoirs
  void s_numFastRes(const unsigned& numRes);//!< set number of fast reservoirs


  void s_parLoadToCalib(std::vector<std::pair<numberSel,par_HRUtype>>& parsToLoad);
  void s_params(const std::pair <numberSel, par_HRUtype>& parDta);
  void s_default();

  void p_param();
  void PDM_boundary_update(); //!< Adjusts the upper and lower parameters cmin and cmax when using the PDM model

  unsigned g_numPars();//!< Get the number of parameters
  void current_param(gs_STORtype gs_STORAGE,soil_STORtype soil_STORAGE,interception_STORtype intrc_STORAGE,surface_STORtype srfs_STORAGE,fast_Response fast_RESP);
  std::vector<std::string> par_HRUtype_to_string(std::list<par_HRUtype> par_list);
  void print_par_list(std::list<par_HRUtype> par_list);


  // std::vector<std::string> Current_parameter_string;
  //std::vector<double> Current_parameter_val;
  //std::vector<double> Current_upparameter_val;
  //std::vector<double> Current_lowparameter_val;

  unsigned g_sizeVecNamesPars();

  std::vector<std::string>  get_CurParNames();
  std::vector<double>  get_CurParVals();
  std::vector<double>  get_CurUpParVals();
  std::vector<double>  get_CurLowParVals();

 protected:

 private:
  //!< The parameters for single HMunit
  unsigned numPars;//!<  The number of model parameters
  pdata pars;//!< The values of parameters on given HRU
  pdata up_pars;//!< The upper bounds of parameters
  pdata low_pars;//!< The lower bound of parameters
  unsigned numFastRes;//!< The number of fast runoff reservoirs, states are implemented in data_HB_1d class

  std::vector<std::string> Current_parameter_string;//!< The vector of param names
  std::list<par_HRUtype> Current_parameter_list;//!<
  std::vector<double> Current_parameter_val;//!< The vector of current parameter values
  std::vector<double> Current_upparameter_val;//!< vector of upper limit values for current parameters
  std::vector<double> Current_lowparameter_val;//!< vector of lower limit values for current parameters

};

#endif // PARAMS_H
