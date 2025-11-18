#ifndef PARSTRUCTSELS_H
#define PARSTRUCTSELS_H

#include <vector>
#include <utility> // pair
#include <string> //string
#include <list>


//srfs_STORAGE
std::list<par_HRUtype> L_SurfaceAll = { par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::TETR, par_HRUtype::RETCAP };
std::list<par_HRUtype> L_SurfacePRTL = {  };
std::list<par_HRUtype> L_Wetland = {  };
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
std::list<par_HRUtype> L_SerialCascadeLinRes = { par_HRUtype::KFR, par_HRUtype::ADIV };
std::list<par_HRUtype> L_SerialLinResGWGros = { par_HRUtype::KFR, par_HRUtype::ADIV,par_HRUtype::RBEI };
std::list<par_HRUtype> L_SerialLinResSoilSois = { par_HRUtype::KFR, par_HRUtype::ADIV,par_HRUtype::RBAI };
std::list<par_HRUtype> L_SerialLinResGWGrosSoilSois = { par_HRUtype::KFR, par_HRUtype::ADIV,par_HRUtype::RBEI, par_HRUtype::RBAI };
//other
std::list<par_HRUtype> L_interception_snow = { par_HRUtype::TETR };
std::list<par_HRUtype> L_snow_melt = { par_HRUtype::DDFA, par_HRUtype::TMEL };


#endif // PARSTRUCTSELS_H
