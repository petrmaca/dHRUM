#ifndef PONDSEL_H_INCLUDED
#define PONDSEL_H_INCLUDED

#include <iomanip>
#include <limits>
#include <string>
#include <map>

//existence of a pond
enum class pond_type {noPond,Pond};
const pond_type all_ponds[] {pond_type::noPond,pond_type::Pond};
const std::vector<std::string> all_pondNames {"noPond","Pond"};

//water surface evaporation variants
enum class ETpond_type {ETpond1, ETpond2};
const ETpond_type all_ETpond[]{ ETpond_type::ETpond1,ETpond_type::ETpond2};
const std::vector<std::string> allETpondTypeNames {"ETpond1", "ETpond2"};

//soil percolation variants
enum class PondSOISPerc_type {noPondSOISPerc, PondSOISPerc1, PondSOISPerc2, PondSOISPerc3};
const PondSOISPerc_type all_PondSOISPerc[]{ PondSOISPerc_type::noPondSOISPerc,PondSOISPerc_type::PondSOISPerc1,PondSOISPerc_type::PondSOISPerc2,PondSOISPerc_type::PondSOISPerc3};
const std::vector<std::string> allPondSOISPondTypeNames {"noPondSOISPerc","PondSOISPerc1", "PondSOISPerc2", "PondSOISPerc3"};

//groundwater percolation variants
enum class PondGWPerc_type {noPondGWPerc,PondGWPerc1, PondGWPerc2,PondGWPerc3};
const PondGWPerc_type all_PondGWPerc[]{ PondGWPerc_type::noPondGWPerc,PondGWPerc_type::PondGWPerc1,PondGWPerc_type::PondGWPerc2,PondGWPerc_type::PondGWPerc3};
const std::vector<std::string> allPondGWPercTypeNames {"noPondGWPerc","PondGWPerc1", "PondGWPerc2","PondGWPerc3"};

//regular out variants
enum class PondRouT_type {noPondRouT,PondRouT1, PondRouT2,PondRouT3};
const PondRouT_type all_PondRouT[]{ PondRouT_type::noPondRouT,PondRouT_type::PondRouT1,PondRouT_type::PondRouT2,PondRouT_type::PondRouT3};
const std::vector<std::string> allPondRouTNames {"noPondRouT","PondRouT1", "PondRouT2","PondRouT3"};

//inputs needed from user
enum class PInp{Area,Max,MRF,CoflW,ET,inSOIS,inGW,outSOIS,outGW,outReg};
const PInp all_pond[]{PInp::Area,PInp::Max,PInp::MRF,PInp::CoflW,PInp::ET,PInp::inSOIS,PInp::inGW,PInp::outSOIS,PInp::outGW,PInp::outReg};
const std::vector<std::string> allPondNames {"Area","Max","MRF","CoflW","ET","inSOIS","inGW","outSOIS","outGW","outReg"};


#endif // PONDSEL_H_INCLUDED
