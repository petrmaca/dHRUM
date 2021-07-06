#ifndef NUMBERSEL_H_INCLUDED
#define NUMBERSEL_H_INCLUDED

#include <iomanip>
#include <limits>
#include <string>

// Type for numbers selector

using numberSel = double;
typedef std::numeric_limits< numberSel > dbl;
using numberDta = unsigned;

using hdata = std::valarray<numberSel>;
using caldata = std::valarray<unsigned>;

enum class ts_type {PREC,SNOW,AET,PET,TEMP,MELT,TROF,STEF,CANF,CANS,STES,EVAC,EVAS,EVBS,INTS,SOIS,GROS,SURS,TOTR,BASF,DIRR,PERC,PREF};
const ts_type all_ts[]{ts_type::PREC,ts_type::SNOW,ts_type::AET,ts_type::PET,ts_type::TEMP,ts_type::MELT,ts_type::TROF,ts_type::STEF,ts_type::CANF, \
                                   ts_type::CANS,ts_type::STES,ts_type::EVAC,ts_type::EVAS,ts_type::EVBS,ts_type::INTS,ts_type::SOIS,ts_type::GROS, \
                                   ts_type::SURS,ts_type::TOTR,ts_type::BASF,ts_type::DIRR,ts_type::PERC,ts_type::PREF};
const unsigned numTSvars = 23;
enum class init_Stype {SOIL, GROUNDWAT, FASTRES, SURFRET, CANS, STES, SNOS};
enum class pet_Type {OUDIN, HAMON};

enum class cal_Type {YEAR, MONTH, DAY, JDAY};
const cal_Type all_caDT[] {cal_Type::YEAR, cal_Type::MONTH, cal_Type::DAY, cal_Type::JDAY};

using pdata = std::valarray<numberSel>;
enum class par_HRUtype {B_SOIL, C_MAX, B_EVAP, SMAX, KS, KF, ADIV, CDIV, SDIV, CAN_ST, STEM_ST, CSDIV, TETR, DDFA, TMEL, RETCAP, L, D_BYPASS, B_EXP};
const par_HRUtype all_pars[]{par_HRUtype::B_SOIL, par_HRUtype::C_MAX, par_HRUtype::B_EVAP, par_HRUtype::SMAX, par_HRUtype::KS, par_HRUtype::KF,      \
                             par_HRUtype::ADIV, par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::CAN_ST, par_HRUtype::STEM_ST, par_HRUtype::CSDIV, \
                             par_HRUtype::TETR, par_HRUtype::DDFA, par_HRUtype::TMEL, par_HRUtype::RETCAP, par_HRUtype::L, par_HRUtype::D_BYPASS, par_HRUtype::B_EXP};

const std::vector<std::string> allParNames {"B_SOIL","C_MAX","B_EVAP","KS","KF","ADIV","CDIV", \
                                      "SDIV","CAN_ST","CAN_ST","STEM_ST","CSDIV","TETR", \
                                      "DDFA","TMEL","RETCAP","L", "D_BYPASS", "B_EXP"};
enum class gs_STORtype { LIN_RES, LINL_RES, LINBY_RES, POW_RES, EXP_RES };
const gs_STORtype all_gs_STORs[]{ gs_STORtype::LIN_RES, gs_STORtype::LINL_RES, gs_STORtype::LINBY_RES, gs_STORtype::POW_RES, gs_STORtype::EXP_RES };

#endif // NUMBERSEL_H_INCLUDED
