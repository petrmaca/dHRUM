#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <valarray>
#include <vector>
#include <utility> // pair
#include <string> //string

#include "numberSel.h"
#include <list>

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

  unsigned g_numPars();//!< Get the number of parameters
  void current_param(gs_STORtype gs_STORAGE,soil_STORtype soil_STORAGE,interception_STORtype intrc_STORAGE,surface_STORtype srfs_STORAGE,fast_Response fast_RESP);

  std::list<par_HRUtype>Current_parameter_list;
  std::vector<double>Current_parameter_val;
  std::vector<double>Current_upparameter_val;
  std::vector<double>Current_lowparameter_val;

 protected:

 private:
  //!< The parameters for single HMunit
  unsigned numPars;//!<  The number of model parameters
  pdata pars;//!< The values of parameters on given HRU
  pdata up_pars;//!< The upper bounds of parameters
  pdata low_pars;//!< The lower bound of parameters
  unsigned numFastRes;//!< The number of fast runoff reservoirs, states are implemented in data_HB_1d class




  //srfs_STORAGE
  std::list<par_HRUtype> L_SurfaceAll = { par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::TETR, par_HRUtype::RETCAP };
  std::list<par_HRUtype> L_SurfacePRTL = {  };
  //intrc_STORAGE
  std::list<par_HRUtype> L_Rutter_Gash = { par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::CAN_ST, par_HRUtype::STEM_ST,par_HRUtype::CSDIV };
  //gs_STORAGE
  std::list<par_HRUtype> L_LIN_RES = { par_HRUtype::KS, par_HRUtype::ADIV };
  std::list<par_HRUtype> L_LINL_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::L };
  std::list<par_HRUtype> L_LINBY_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::D_BYPASS };
  std::list<par_HRUtype> L_POW_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::B_EXP };
  std::list<par_HRUtype> L_EXP_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::B_EXP };
  std::list<par_HRUtype> L_LIN_2SE = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::KS2 };
  std::list<par_HRUtype> L_LIN_2PA = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::KS2, par_HRUtype::ALPHA };
  std::list<par_HRUtype> L_FLEX_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::KS2, par_HRUtype::THR };
  std::list<par_HRUtype> L_EXP_LOG = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::B_EXP };
  //soil_STORAGE
  std::list<par_HRUtype> L_PDM = { par_HRUtype::B_SOIL, par_HRUtype::C_MAX, par_HRUtype::B_EVAP, par_HRUtype::SMAXpdm, par_HRUtype::CMIN };
  std::list<par_HRUtype> L_COLLIE_V2 = { par_HRUtype::KF, par_HRUtype::FC, par_HRUtype::FOREST_FRACT, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_NEW_ZEALAND = { par_HRUtype::KF, par_HRUtype::FC, par_HRUtype::FOREST_FRACT, par_HRUtype::KF2, par_HRUtype::KF_NONLIN, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_GR4J = { par_HRUtype::SMAX};
  std::list<par_HRUtype> L_SBROOK_V1 = { par_HRUtype::KF, par_HRUtype::FC, par_HRUtype::FOREST_FRACT, par_HRUtype::KF_NONLIN, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_HILLSLOPE = { par_HRUtype::KF, par_HRUtype::KF_NONLIN, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_PLATEAU = { par_HRUtype::C, par_HRUtype::INFR_MAX, par_HRUtype::RF, par_HRUtype::WP, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_PDM2 = { par_HRUtype::B_SOIL, par_HRUtype::C_MAX };
  //fast_response
  std::list<par_HRUtype> L_SerialCascadeLinRes = { par_HRUtype::KF, par_HRUtype::ADIV };
  std::list<par_HRUtype> L_SerialLinResGWGros = { par_HRUtype::KF, par_HRUtype::ADIV,par_HRUtype::RBEI };
  std::list<par_HRUtype> L_SerialLinResSoilSois = { par_HRUtype::KF, par_HRUtype::ADIV,par_HRUtype::RBAI };
  std::list<par_HRUtype> L_SerialLinResGWGrosSoilSois = { par_HRUtype::KF, par_HRUtype::ADIV,par_HRUtype::RBEI, par_HRUtype::RBAI };
  //other
  std::list<par_HRUtype> L_interception_snow = { par_HRUtype::TETR };
  std::list<par_HRUtype> L_snow_melt = { par_HRUtype::DDFA, par_HRUtype::TMEL };
};

#endif // PARAMS_H
