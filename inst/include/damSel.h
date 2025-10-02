#ifndef DAMSEL_H
#define DAMSEL_H

#include <iomanip>
#include <limits>
#include <string>
#include <map>


const unsigned damTSvars = 10;
enum class dam_ts {PREC,DAMS,INFL,DAIS,DAIG,OUFL,ETDM,OFLW,OULT,INLT};
const dam_ts all_dts[]{dam_ts::PREC,dam_ts::DAMS,dam_ts::INFL,dam_ts::DAIS,dam_ts::DAIG,\
                       dam_ts::OUFL,dam_ts::ETDM,dam_ts::OFLW,dam_ts::OULT,dam_ts::INLT};

enum class init_dStype {DAMS};



//water surface evaporation variants
enum class ETdam_type {ETdam1, ETdam2};
const ETdam_type all_ETdam[]{ ETdam_type::ETdam1,ETdam_type::ETdam2};
const std::vector<std::string> allETdamTypeNames {"ETpond1, ETpond2"};

//soil percolation variants
enum class DamSOISPerc_type {noDamSOISPerc, DamSOISPerc1, DamSOISPerc2, DamSOISPerc3};
const DamSOISPerc_type all_DamSOISPerc[]{ DamSOISPerc_type::noDamSOISPerc,DamSOISPerc_type::DamSOISPerc1,DamSOISPerc_type::DamSOISPerc2,DamSOISPerc_type::DamSOISPerc3};
const std::vector<std::string> allDamSOISPercTypeNames {"noDamSOISPerc","DamSOISPerc1", "DamSOISPerc2", "DamSOISPerc3"};

//groundwater percolation variants
enum class DamGWPerc_type {noDamGWPerc,DamGWPerc1, DamGWPerc2,DamGWPerc3};
const DamGWPerc_type all_DamGWPerc[]{ DamGWPerc_type::noDamGWPerc,DamGWPerc_type::DamGWPerc1,DamGWPerc_type::DamGWPerc2,DamGWPerc_type::DamGWPerc3};
const std::vector<std::string> allDamGWPercTypeNames {"noDamGWPerc","DamGWPerc1", "DamGWPerc2","DamGWPerc3"};

//regular out variants
enum class DamRouT_type {noDamRouT,DamRouT1, DamRouT2,DamRouT3};
const DamRouT_type all_DamRouT[]{ DamRouT_type::noDamRouT,DamRouT_type::DamRouT1,DamRouT_type::DamRouT2,DamRouT_type::DamRouT3};
const std::vector<std::string> allDamRouTNames {"noDamRouT","DamRouT1", "DamRouT2","DamRouT3"};



#endif // DAMSEL_H
