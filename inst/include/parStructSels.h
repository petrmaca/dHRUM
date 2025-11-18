#ifndef PARSTRUCTSELS_H
#define PARSTRUCTSELS_H

#include <vector>
#include <utility> // pair
#include <string> //string
#include <list>

#include "numberSel.h"

enum class pet_Type {OUDIN, HAMON, THORNTHWAITE, BLANEYCRIDDLE, JENSENHAISE, MCGUINNESSBORDNE};

const std::vector<std::string> allPetNames {"OUDIN", "HAMON", "THORNTHWAITE","BLANEYCRIDDLE","JENSENHAISE", "MCGUINNESSBORDNE"};

using pdata = std::valarray<numberSel>;
enum class par_HRUtype {B_SOIL, C_MAX, B_EVAP, SMAXpdm, KS, KF, ADIV, CDIV, SDIV, CAN_ST, STEM_ST, CSDIV, TETR, DDFA, TMEL, RETCAP, L, D_BYPASS, B_EXP, KS2, THR, ALPHA, CMIN, FC, FOREST_FRACT, KF2, KF_NONLIN, C, INFR_MAX, RF, WP, SMAX, RBAI, RBEI,KFR};
const par_HRUtype all_pars[]{par_HRUtype::B_SOIL, par_HRUtype::C_MAX, par_HRUtype::B_EVAP, par_HRUtype::SMAXpdm, par_HRUtype::KS, par_HRUtype::KF,      \
                             par_HRUtype::ADIV, par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::CAN_ST, par_HRUtype::STEM_ST, par_HRUtype::CSDIV,    \
                             par_HRUtype::TETR, par_HRUtype::DDFA, par_HRUtype::TMEL, par_HRUtype::RETCAP, par_HRUtype::L, par_HRUtype::D_BYPASS,
                             par_HRUtype::B_EXP, par_HRUtype::KS2, par_HRUtype::THR, par_HRUtype::ALPHA, par_HRUtype::CMIN, \
                             par_HRUtype::FC, par_HRUtype::FOREST_FRACT, par_HRUtype::KF2, par_HRUtype::KF_NONLIN,          \
                             par_HRUtype::C, par_HRUtype::INFR_MAX, par_HRUtype::RF, par_HRUtype::WP, par_HRUtype::SMAX,    \
                             par_HRUtype::RBAI,par_HRUtype::RBEI,par_HRUtype::KFR};

const std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP", "SMAXpdm","KS","KF","ADIV","CDIV", \
                                            "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR",            \
                                            "DDFA","TMEL","RETCAP","L", "D_BYPASS", "B_EXP", "KS2", "THR",
                                            "ALPHA","CMIN","FC","FOREST_FRACT", "KF2", "KF_NONLIN", "C", \
                                            "INFR_MAX", "RF", "WP", "SMAX", "RBAI", "RBEI", "KFR"};

enum class gs_STORtype { LIN_RES, LINL_RES, LINBY_RES, POW_RES, EXP_RES, LIN_2SE, LIN_2PA, FLEX_RES,EXP_LOG};
const gs_STORtype all_gs_STORs[]{ gs_STORtype::LIN_RES, gs_STORtype::LINL_RES, gs_STORtype::LINBY_RES,
                                  gs_STORtype::POW_RES, gs_STORtype::EXP_RES, gs_STORtype::LIN_2SE,
                                  gs_STORtype::LIN_2PA, gs_STORtype::FLEX_RES, gs_STORtype::EXP_LOG };

const std::vector<std::string> allGWStorTypeNames {"LIN_RES", "LINL_RES", "LINBY_RES", "POW_RES",
                                                   "EXP_RES", "LIN_2SE", "LIN_2PA", "FLEX_RES","EXP_LOG" };

enum class soil_STORtype { PDM, COLLIE_V2, NEW_ZEALAND, GR4J, SBROOK_V1, HILLSLOPE, PLATEAU, PDM2 };
const soil_STORtype all_soil_STORs[]{ soil_STORtype::PDM, soil_STORtype::COLLIE_V2, soil_STORtype::NEW_ZEALAND, \
                                      soil_STORtype::GR4J, soil_STORtype::SBROOK_V1, soil_STORtype::HILLSLOPE,  \
                                      soil_STORtype::PLATEAU,soil_STORtype::PDM2};

const std::vector<std::string> allSoilStorTypeNames {"PDM", "COLLIE_V2", "NEW_ZEALAND", "GR4J", "SBROOK_V1", \
                                                     "HILLSLOPE", "PLATEAU","PDM2"};

enum class interception_STORtype { Rutter_Gash};
const interception_STORtype all_Interceptions[]{ interception_STORtype::Rutter_Gash};
const std::vector<std::string> allTnterceptionStorTypeNames {"Rutter_Gash"};

enum class surface_STORtype {SurfaceAll, SurfacePRTL, Wetland};
const surface_STORtype all_SurfacStors[] {surface_STORtype::SurfaceAll, surface_STORtype::SurfacePRTL,surface_STORtype::Wetland};
const std::vector<std::string> allSurfacStorsNames {"SurfaceAll","SurfacePRTL","Wetland"};

enum class fast_Response {SerialCascadeLinRes, SerialLinResGWGros,SerialLinResSoilSois,SerialLinResGWGrosSoilSois};
const fast_Response all_FastResponses[] {fast_Response::SerialCascadeLinRes, fast_Response::SerialLinResGWGros, fast_Response::SerialLinResSoilSois, fast_Response::SerialLinResGWGrosSoilSois};
const std::vector<std::string> all_FastResponsesNames {"SerialCascadeLinRes", "SerialLinResGWGros","SerialLinResSoilSois", "SerialLinResGWGrosSoilSois"};


//srfs_STORAGE
const std::list<par_HRUtype> L_SurfaceAll = { par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::TETR, par_HRUtype::RETCAP };
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
