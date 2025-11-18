#ifndef NUMBERSEL_H_INCLUDED
#define NUMBERSEL_H_INCLUDED

#include <iomanip>
#include <limits>
#include <string>
#include <map>

#include "parStructSels.h"

// Type for numbers selector

using numberSel = double;
typedef std::numeric_limits< numberSel > dbl;
using numberDta = unsigned;

using hdata = std::valarray<numberSel>;
using caldata = std::valarray<unsigned>;

enum class ts_type {PREC,SNOW,AET,PET,TEMP,MELT,TROF,STEF,CANF,CANS,STES,EVAC,EVAS,EVBS,INTS,SOIS,GROS,SURS,TOTR,BASF,DIRR,PERC,PREF,ETSW,PONS,ETPO,POIS,POIG};
const ts_type all_ts[]{ts_type::PREC,ts_type::SNOW,ts_type::AET,ts_type::PET,ts_type::TEMP,ts_type::MELT,ts_type::TROF,ts_type::STEF,ts_type::CANF, \
                                   ts_type::CANS,ts_type::STES,ts_type::EVAC,ts_type::EVAS,ts_type::EVBS,ts_type::INTS,ts_type::SOIS,ts_type::GROS, \
                                   ts_type::SURS,ts_type::TOTR,ts_type::BASF,ts_type::DIRR,ts_type::PERC,ts_type::PREF,ts_type::ETSW,ts_type::PONS,\
                                   ts_type::ETPO,ts_type::POIS,ts_type::POIG};
//ts_type::ETPO,ts_type::POIS,ts_type::POIG

const unsigned numTSvars = 28;
enum class init_Stype {SOIL, GROUNDWAT, FASTRES, SURFRET, CANS, STES, SNOS,GROS1,GROS2};

enum class cal_Type {YEAR, MONTH, DAY, JDAY};
const cal_Type all_caDT[] {cal_Type::YEAR, cal_Type::MONTH, cal_Type::DAY, cal_Type::JDAY};



#endif // NUMBERSEL_H_INCLUDED
