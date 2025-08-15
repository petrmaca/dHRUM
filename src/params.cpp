#include "params.h"

params::params(): numPars(0),
  pars(1,1),
  up_pars(1,1),
  low_pars(1,1),
  numFastRes(1) {
  //ctor
//    b_soil = 2.0;
//    c_max = 100.0;
//    b_evap = 1;
  numPars = 34;

  pars.resize(numPars,numPars);
  up_pars.resize(numPars,numPars);
  low_pars.resize(numPars,numPars);
  pars[0] = 2.0;//!< B_SOIL Parameter controlling shape of Pareto distribution of soil storages [0,inf] however [0.5,3],VC1
  pars[1] = 10.0;//!< C_MAX Max storage of storages distributed by Pareto distribution [0,inf],VC1
  pars[2] = 1.0;//!< B_EVAP Parameter controlling soil evapotranspiration [0,infty] how ever [0.5,3],VC1
  pars[3] = 0.0001;//!< SMAXpdm Max soil storage calculate using Cmax and b_soil
  pars[4] = 0.1;//!< KS Storage coefficient of groundwater storage [0,1],VC1
  pars[5] = 0.5;//!< KF Storage coefficient of runoff response reservoirs [0,1],VC1
  pars[6] = 0.5;//!< ADIV Divider of percolation into the direct runoff input par[6]*Perc and  groundwater input (1-par[6])*Perc [0,1],VC1
  pars[7] = 0.02;//!< CDIV Divider of gross rainfall as a Canopy input [0,1],VC1
  pars[8] = 0.02;//!< SDIV Divider of gross rainfall as a Trunk input [0,1],VC1
  pars[9] = 1.0;//!< CAN_ST The Max canopy storage [0,inf],VC1
  pars[10] = 1.0;//!< STEM_ST The Max stem and trunk storage [0,inf],VC1
  pars[11] = 0.5;//!< CSDIV The divider of canopy outflow to throughflow and stemflow storage [0,1],VC1
  pars[12] = 1;//!< TETR The threshold temperature for determining snow [-inf,inf] better [-5,5]
  pars[13] = 2;//!<  DDFA The day degree model for snow melt [o, inf] better [0,2],VC1
  pars[14] = 0;//!<  TMEL The threshold temperature for determining melting process [-inf, inf] better [-5,5]
  pars[15] = 4;//!< RETCAP The maximum capacity of surface retention [0, inf],VC1
  pars[16] = 1;//!< L The amount of groundwater recharge removed from the linear reservoir [0,1]
  pars[17] = 1;//!< D_BYPASS The amount of groundwater recharge removed from the linear reservoir [0,1]
  pars[18] = 1;//!< B_EXP Power coefficient
  pars[19] = 0.1;//!< KS2 Storage coefficient of groundwater storage [0,1],VC1
  pars[20] = 1;//!< THR Threshold coefficient for threshold-controlled linear storage [0,inf]
  pars[21] = 0.5;//!< ALPHA Divider for two parallel linear reservoirs
  pars[22] = 0;//!< Cmin for pdmsoil reservoir [0,inf]
  numberSel helpSmax;
  helpSmax = (pars[0] * pars[22] +pars[1]) / (pars[0]+1);
  pars[23] = 10;//!< FC Field capacity [mm] [0,inf)
  pars[24] = 0.5;//!< Forest fraction [0,1]
  pars[25] = 0.5;//!< KF2 Storage coefficient of runoff response reservoirs [0,1]
  pars[26] = 0.5;//!< KF_NONLIN runoff non-linearity parameter [-] [0,inf)
  pars[27] = 1;//!< C capillary rise [mm d-1] [0, inf)
  pars[28] = 1;//!< INFR_MAX Maximum infiltration rate [mm/d], [0,inf)
  pars[29] = 0.5;//!< RF evaporation reduction factor [-] [0,1]
  pars[30] = 0.5;//!< WP wilting point [-] [0,1]
  pars[31] = 20;//!< SMAX  [mm] [0,inf]
  pars[32] = 0.05;//!< RBAI River bank infiltration rate (infiltration to Soil storage)
  pars[33] = 0.1;//!< RBEI River bed infiltration rate (infiltration to ground water storage)
// Upper bounds of parameters
  up_pars[0] = 3.0;//!< B_SOIL Parameter controlling shape of Pareto distribution of soil storages [0,inf] however [0.5,3],VC1
  up_pars[1] = 500.0;//!< C_MAX Max storage of storages distributed by Pareto distribution [0,inf],VC1
  up_pars[2] = 3.0;//!< B_EVAP Parameter controlling soil evapotranspiration [0,infty] how ever [0.5,3],VC1
  up_pars[3] = 0.0;//!< SMAXpdm Max soil storage calculate using Cmax and b_soil
  up_pars[4] = 1;//!< KS Storage coefficient of groundwater storage [0,1],VC1
  up_pars[5] = 1;//!< KF Storage coefficient of runoff response reservoirs [0,1],VC1
  up_pars[6] = 1;//!< ADIV Divider of percolation into the overflow input par[6]*Perc and  groundwater input (1-par[6])*Perc [0,1],VC1
  up_pars[7] = 1;//!< CDIV Divider of gross rainfall as a Canopy input [0,1],VC1
  up_pars[8] = 1;//!< SDIV Divider of gross rainfall as a Trunk input [0,1],VC1
  up_pars[9] = 100;//!< CAN_ST The Max canopy storage [0,inf],VC1
  up_pars[10] = 100;//!< STEM_ST The Max stem and trunk storage [0,inf],VC1
  up_pars[11] = 1;//!< CSDIV The divider of canopy outflow to throughflow and stemflow storage [0,1],VC1
  up_pars[12] = 1;//!< TETR The threshold temperature for determining snow [-inf,inf] better [-5,5]
  up_pars[13] = 2;//!<  DDFA The day degree model for snow melt [o, inf] better [0,2],VC1
  up_pars[14] = 0;//!<  TMEL The threshold temperature for determining melting process [-inf, inf] better [-5,5]
  up_pars[15] = 100;//!< RETCAP The maximum capacity of surface retention [0, inf],VC1
  up_pars[16] = 1;//!< L The amount of groundwater recharge removed from the linear reservoir [0,1]
  up_pars[17] = 1;//!< D_BYPASS The amount of groundwater recharge removed from the linear reservoir [0,1]
  up_pars[18] = 1;//!< B_EXP Power coefficient
  up_pars[19] = 1;//!< KS2 Storage coefficient of groundwater storage [0,1],VC1
  up_pars[20] = 100;//!< THR Threshold coefficient for threshold-controlled linear storage [1,inf]
  up_pars[21] = 1;//!< ALPHA Divider for two parallel linear reservoirs
  up_pars[22] = 200;//!< CIMN lower limit of c in soils pdm reservoir

  numberSel helpSmaxPDMUp;
  helpSmaxPDMUp = (up_pars[0]*up_pars[22] + up_pars[1]) / (up_pars[0]+1);
  pars[3] = helpSmaxPDMUp;

  up_pars[23] = 100;//!< FC Field capacity [mm] [0,inf)
  up_pars[24] = 1;//!< Forest fraction [0,1]
  up_pars[25] = 1;//!< KF2 Storage coefficient of runoff response reservoirs [0,1],VC1
  up_pars[26] = 1;//!< KF_NONLIN runoff non-linearity parameter [-] [0,inf)
  up_pars[27] = 100;//!< C capillary rise [mm d-1] [0, inf)
  up_pars[28] = 100;//!< INFR_MAX Maximum infiltration rate [mm/d], [0,inf)
  up_pars[29] = 1;//!< RF evaporation reduction factor [-] [0,1]
  up_pars[30] = 1;//!< WP wilting point [-] [0,1]
  up_pars[31] = 100;//!< SMAX  [mm] [0,inf]
  up_pars[32] = 1;//!< RBAI River bank infiltration rate (infiltration to Soil storage)
  up_pars[33] = 1;//!< RBEI River bed infiltration rate (infiltration to ground water storage)

// Lower bounds of parameters
  low_pars[0] = 0.0;//!< B_SOIL Parameter controlling shape of Pareto distribution of soil storages [0,inf] however [0.5,3],VC1
  low_pars[1] = 0.0;//!< C_MAX Max storage of storages distributed by Pareto distribution [0,inf],VC1
  low_pars[2] = 0.0;//!< B_EVAP Parameter controlling soil evapotranspiration [0,infty] how ever [0.5,3],VC1
  low_pars[3] = 0.0;//!< SMAXpdm Max soil storage calculate using Cmax and b_soil
  low_pars[4] = 0.0;//!< KS Storage coefficient of groundwater storage [0,1],VC1
  low_pars[5] = 0.0;//!< KF Storage coefficient of runoff response reservoirs [0,1],VC1
  low_pars[6] = 0.0;//!< ADIV Divider of percolation into the overflow input par[6]*Perc and  groundwater input (1-par[6])*Perc [0,1],VC1
  low_pars[7] = 0.0;//!< CDIV Divider of gross rainfall as a Canopy input [0,1],VC1
  low_pars[8] = 0.0;//!< SDIV Divider of gross rainfall as a Trunk input [0,1],VC1
  low_pars[9] = 0.0;//!< CAN_ST The Max canopy storage [0,inf],VC1
  low_pars[10] = 0.0;//!< STEM_ST The Max stem and trunk storage [0,inf],VC1
  low_pars[11] = 0.0;//!< CSDIV The divider of canopy outflow to throughflow and stemflow storage [0,1],VC1
  low_pars[12] = 0.0;//!< TETR The threshold temperature for determining snow [-inf,inf] better [-5,5]
  low_pars[13] = 0.0;//!<  DDFA The day degree model for snow melt [o, inf] better [0,2],VC1
  low_pars[14] = 0.0;//!<  TMEL The threshold temperature for determining melting process [-inf, inf] better [-5,5]
  low_pars[15] = 0.0;//!< RETCAP The maximum capacity of surface retention [0, inf],VC1
  low_pars[16] = 0.0;//!< L The amount of groundwater recharge removed from the linear reservoir [0,1]
  low_pars[17] = 0.0;//!< D_BYPASS The amount of groundwater recharge removed from the linear reservoir [0,1]
  low_pars[18] = 0.25;//!< B_EXP Power coefficient
  low_pars[19] = 0.0;//!< KS2 Storage coefficient of groundwater storage [0,1],VC1
  low_pars[20] = 0.0;//!< THR Threshold coefficient for threshold-controlled linear storage [1,inf]
  low_pars[21] = 0.0;//!< ALPHA Divider for two parallel linear reservoirs
  low_pars[22] = 0;//!< CMIN lower limit of c in soils pdm reservoir

  numberSel helpSmaxPDMLow;
  helpSmaxPDMLow = (low_pars[0] *low_pars[22] +low_pars[1]) / (low_pars[0]+1);
  low_pars[3] = helpSmaxPDMLow;

  low_pars[23] = 0;//!< FC Field capacity [mm] [0,inf)
  low_pars[24] = 0.0;//!< Forest fraction [0,1]
  low_pars[25] = 0.0;//!< KF2 Storage coefficient of runoff response reservoirs [0,1]
  low_pars[26] = 0.0;//!< KF_NONLIN runoff non-linearity parameter [-] [0,inf)
  low_pars[27] = 0.0;//!< C capillary rise [mm d-1] [0, inf)
  low_pars[28] = 0.0;//!< INFR_MAX Maximum infiltration rate [mm/d], [0,inf)
  low_pars[29] = 0.0;//!< RF evaporation reduction factor [-] [0,1]
  low_pars[30] = 0.0;//!< WP wilting point [-] [0,1]
  low_pars[31] = 0.0;//!< SMAX  [mm] [0,inf]
  low_pars[32] = 0;//!< RBAI River bank infiltration rate (infiltration to Soil storage)
  low_pars[33] = 0;//!< RBEI River bed infiltration rate (infiltration to ground water storage)

//  std::cout << "Params are initialized." << std::endl;
}

params::~params() {
  //dtor
}

params::params(const params& other): numPars(0),
  pars(1,1),
  up_pars(1,1),
  low_pars(1,1),
  numFastRes(1)  {

  numPars = other.numPars;
  pars = other.pars;
  up_pars = other.up_pars;
  low_pars = other.low_pars;
  numFastRes = other.numFastRes;

}

params& params::operator=(const params& rhs) {

  if (this == &rhs) return *this;
  else {
    numPars = rhs.numPars;
    pars = rhs.pars;
    up_pars = rhs.up_pars;
    low_pars = rhs.low_pars;
    numFastRes = rhs.numFastRes;
  }

  return *this;
}

/** \brief Setting the value of parameter using the pair, overloaded function
 *
 * \param value of param contorled by its par_HRUtype
 * \param ID of selected parameter controlled by the par_HRUtype
 *
 */
void params::s_params(const numberSel& par_dta,par_HRUtype _parType) {

  switch(_parType) {
  case par_HRUtype::B_SOIL:
    pars[0] = par_dta;
    pars[3] = (pars[0]*pars[22] + pars[1]) / (pars[0]+1);
//    std::cout << "New b_soil --> loaded\n";
    break;
  case par_HRUtype::C_MAX:
    pars[1] = par_dta;
    pars[3] = (pars[0]*pars[22] +  pars[1]) / (pars[0]+1);
//    std::cout << "New c_max --> loaded\n";
    break;
  case par_HRUtype::B_EVAP:
    pars[2] = par_dta;
//    std::cout << "New b_evap --> loaded\n";
    break;
  case par_HRUtype::SMAXpdm:
    pars[3] = par_dta;
    pars[1] = pars[3]*(pars[0]+1) - (pars[0]*pars[22]);//For preventing the losing the link between b_soil and cmax and cmin see pdm paper hess 2007
//    std::cout << "New Smax --> loaded\n";
    break;
  case par_HRUtype::KS:
    pars[4] = par_dta;
//    std::cout << "New Ks --> loaded\n";
    break;
  case par_HRUtype::KF:
    pars[5] = par_dta;
//    std::cout << "New Kf --> loaded\n";
    break;
  case par_HRUtype::ADIV:
    pars[6] = par_dta;
//    std::cout << "New Adiv --> loaded\n";
    break;
  case par_HRUtype::CDIV:
    pars[7] = par_dta;
//    std::cout << "New Cdiv --> loaded\n";
    break;
  case par_HRUtype::SDIV:
    pars[8] = par_dta;
//    std::cout << "New Sdiv --> loaded\n";
    break;
  case par_HRUtype::CAN_ST:
    pars[9] = par_dta;
//    std::cout << "New Can_St --> loaded\n";
    break;
  case par_HRUtype::STEM_ST:
    pars[10] = par_dta;
//    std::cout << "New Stem_St --> loaded\n";
    break;
  case par_HRUtype::CSDIV:
    pars[11] = par_dta;
//    std::cout << "New CSdiv --> loaded\n";
    break;
  case par_HRUtype::TETR:
    pars[12] = par_dta;
//    std::cout << "New TETR --> loaded\n";
    break;
  case par_HRUtype::DDFA:
    pars[13] = par_dta;
//    std::cout << "New DDFA --> loaded\n";
    break;
  case par_HRUtype::TMEL:
    pars[14] = par_dta;
//    std::cout << "New TETR --> loaded\n";
    break;
  case par_HRUtype::RETCAP:
    pars[15] = par_dta;
//    std::cout << "New RETCAP --> loaded\n";
    break;
  case par_HRUtype::L:
    pars[16] = par_dta;
//    std::cout << "New L --> loaded\n";
    break;
  case par_HRUtype::D_BYPASS:
    pars[17] = par_dta;
//    std::cout << "New D_BYPASS --> loaded\n";
    break;
  case par_HRUtype::B_EXP:
    pars[18] = par_dta;
//    std::cout << "New B_EXP --> loaded\n";
    break;
  case par_HRUtype::KS2:
    pars[19] = par_dta;
//    std::cout << "New KS2 --> loaded\n";
    break;
  case par_HRUtype::THR:
    pars[20] = par_dta;
//    std::cout << "New THR --> loaded\n";
    break;
  case par_HRUtype::ALPHA:
    pars[21] = par_dta;
    break;
  case par_HRUtype::CMIN:
    pars[22] = par_dta;
    pars[3] = (pars[0] * pars[22] +pars[2]) / (pars[0] +1 );//for preventing consisntency beteewn cmin cmax smax in pdm model
//    std::cout << "New ALPHA --> loaded\n";
    break;
  case par_HRUtype::FC:
    pars[23] = par_dta;
    //    std::cout << "New FC --> loaded\n";
    break;
  case par_HRUtype::FOREST_FRACT:
    pars[24] = par_dta;
    //    std::cout << "New FOREST_FRACT --> loaded\n";
    break;
  case par_HRUtype::KF2:
    pars[25] = par_dta;
    //    std::cout << "New KF2 --> loaded\n";
    break;
  case par_HRUtype::KF_NONLIN:
    pars[26] = par_dta;
    //    std::cout << "New KF_NONLIN --> loaded\n";
    break;
  case par_HRUtype::C:
    pars[27] = par_dta;
    //    std::cout << "New C --> loaded\n";
    break;
  case par_HRUtype::INFR_MAX:
    pars[28] = par_dta;
    //    std::cout << "New INFR_MAX --> loaded\n";
    break;
  case par_HRUtype::RF:
    pars[29] = par_dta;
    //    std::cout << "New RF --> loaded\n";
    break;
  case par_HRUtype::WP:
    pars[30] = par_dta;
    //    std::cout << "New WP --> loaded\n";
    break;
  case par_HRUtype::SMAX:
    pars[31] = par_dta;
    break;
  case par_HRUtype::RBAI:
    pars[32] = par_dta;
    break;
  case par_HRUtype::RBEI:
    pars[33] = par_dta;
    break;
  }

  pars[3] = (pars[0] * pars[22] +pars[2]) / (pars[0] +1 );

  return ;
}

/** \brief Setting the value of parameter using the pair, overloaded function
 *
 * \param pair of value of param and it ID controlled by the par_HRUtype
 *
 */
void params::s_params(const std::pair <numberSel,par_HRUtype>& parDta) {

  numberSel par_dta = parDta.first;
  par_HRUtype _parType = parDta.second;

  switch(_parType) {
  case par_HRUtype::B_SOIL:
    pars[0] = par_dta;
    pars[3] = (pars[1] + pars[0] * pars[22])/ (pars[0]+1);
//    std::cout << "New b_soil --> loaded\n";
    break;
  case par_HRUtype::C_MAX:
    pars[1] = par_dta;
    pars[3] = (pars[1] + pars[0] * pars[22])/ (pars[0]+1);
//    std::cout << "New c_max --> loaded\n";
    break;
  case par_HRUtype::B_EVAP:
    pars[2] = par_dta;
//    std::cout << "New b_evap --> loaded\n";
    break;
  case par_HRUtype::SMAXpdm:
    pars[3] = par_dta;
    pars[1] = pars[3]*(pars[0]+1) -(pars[0]*pars[22]);
//    std::cout << "New Smax --> loaded\n";
    break;
  case par_HRUtype::KS:
    pars[4] = par_dta;
//    std::cout << "New Ks --> loaded\n";
    break;
  case par_HRUtype::KF:
    pars[5] = par_dta;
//    std::cout << "New Kf --> loaded\n";
    break;
  case par_HRUtype::ADIV:
    pars[6] = par_dta;
//    std::cout << "New Adiv --> loaded\n";
    break;
  case par_HRUtype::CDIV:
    pars[7] = par_dta;
//    std::cout << "New Cdiv --> loaded\n";
    break;
  case par_HRUtype::SDIV:
    pars[8] = par_dta;
//    std::cout << "New Sdiv --> loaded\n";
    break;
  case par_HRUtype::CAN_ST:
    pars[9] = par_dta;
//    std::cout << "New Can_St --> loaded\n";
    break;
  case par_HRUtype::STEM_ST:
    pars[10] = par_dta;
//    std::cout << "New Stem_St --> loaded\n";
    break;
  case par_HRUtype::CSDIV:
    pars[11] = par_dta;
//    std::cout << "New CSdiv --> loaded\n";
    break;
  case par_HRUtype::TETR:
    pars[12] = par_dta;
//    std::cout << "New TETR --> loaded\n";
    break;
  case par_HRUtype::DDFA:
    pars[13] = par_dta;
//    std::cout << "New DDFA --> loaded\n";
    break;
  case par_HRUtype::TMEL:
    pars[14] = par_dta;
//    std::cout << "New TETR --> loaded\n";
    break;
  case par_HRUtype::RETCAP:
    pars[15] = par_dta;
//    std::cout << "New RETCAP --> loaded\n";
    break;
  case par_HRUtype::L:
    pars[16] = par_dta;
//    std::cout << "New L --> loaded\n";
    break;
  case par_HRUtype::D_BYPASS:
    pars[17] = par_dta;
//    std::cout << "New D_BYPASS --> loaded\n";
    break;
  case par_HRUtype::B_EXP:
    pars[18] = par_dta;
//    std::cout << "New B_EXP --> loaded\n";
    break;
  case par_HRUtype::KS2:
    pars[19] = par_dta;
//    std::cout << "New KS2 --> loaded\n";
    break;
  case par_HRUtype::THR:
    pars[20] = par_dta;
//    std::cout << "New THR --> loaded\n";
    break;
  case par_HRUtype::ALPHA:
    pars[21] = par_dta;
//    std::cout << "New ALPHA --> loaded\n";
    break;
  case par_HRUtype::CMIN:
    pars[22] = par_dta;
    pars[3] = (pars[1] + pars[0] * pars[22])/ (pars[0]+1);
//    std::cout << "New ALPHA --> loaded\n";
    break;
  case par_HRUtype::FC:
    pars[23] = par_dta;
    //    std::cout << "New FC --> loaded\n";
    break;
  case par_HRUtype::FOREST_FRACT:
    pars[24] = par_dta;
    //    std::cout << "New FOREST_FRACT --> loaded\n";
    break;
  case par_HRUtype::KF2:
    pars[25] = par_dta;
    //    std::cout << "New KF2 --> loaded\n";
    break;
  case par_HRUtype::KF_NONLIN:
    pars[26] = par_dta;
    //    std::cout << "New KF_NONLIN --> loaded\n";
    break;
  case par_HRUtype::C:
    pars[27] = par_dta;
    //    std::cout << "New C --> loaded\n";
    break;
  case par_HRUtype::INFR_MAX:
    pars[28] = par_dta;
    //    std::cout << "New INFR_MAX --> loaded\n";
    break;
  case par_HRUtype::RF:
    pars[29] = par_dta;
    //    std::cout << "New RF --> loaded\n";
    break;
  case par_HRUtype::WP:
    pars[30] = par_dta;
    //    std::cout << "New WP --> loaded\n";
    break;
  case par_HRUtype::SMAX:
    pars[31] = par_dta;
    break;
  case par_HRUtype::RBAI:
    pars[32] = par_dta;
    break;
  case par_HRUtype::RBEI:
    pars[33] = par_dta;
    break;
  }
  return ;
}

/** \brief Getter for param value
 *
 * \param par ID given by par_HRUtype
 *
 * \return value of par controled by its ID par_HRUtype
 *
 */
numberSel params::g_par(const par_HRUtype& _parType) {

  numberSel value = 0.0;

  switch(_parType) {
  case par_HRUtype::B_SOIL:
    value = pars[0];
    break;
  case par_HRUtype::C_MAX:
    value =  pars[1];
    break;
  case par_HRUtype::B_EVAP:
    value =  pars[2];
    break;
  case par_HRUtype::SMAXpdm:
    value =  pars[3];
    break;
  case par_HRUtype::KS:
    value =  pars[4];
    break;
  case par_HRUtype::KF:
    value =  pars[5];
    break;
  case par_HRUtype::ADIV:
    value =  pars[6];
    break;
  case par_HRUtype::CDIV:
    value =  pars[7];
    break;
  case par_HRUtype::SDIV:
    value =  pars[8];
    break;
  case par_HRUtype::CAN_ST:
    value =  pars[9];
    break;
  case par_HRUtype::STEM_ST:
    value =  pars[10];
    break;
  case par_HRUtype::CSDIV:
    value =  pars[11];
    break;
  case par_HRUtype::TETR:
    value =  pars[12];
    break;
  case par_HRUtype::DDFA:
    value =  pars[13];
    break;
  case par_HRUtype::TMEL:
    value =  pars[14];
    break;
  case par_HRUtype::RETCAP:
    value =  pars[15];
    break;
  case par_HRUtype::L:
    value =  pars[16];
    break;
  case par_HRUtype::D_BYPASS:
    value =  pars[17];
    break;
  case par_HRUtype::B_EXP:
    value =  pars[18];
    break;
  case par_HRUtype::KS2:
    value =  pars[19];
    break;
  case par_HRUtype::THR:
    value =  pars[20];
    break;
  case par_HRUtype::ALPHA:
    value =  pars[21];
    break;
  case par_HRUtype::CMIN:
    value =  pars[22];
    break;
  case par_HRUtype::FC:
    value =  pars[23];
    break;
  case par_HRUtype::FOREST_FRACT:
    value =  pars[24];
    break;
  case par_HRUtype::KF2:
    value =  pars[25];
    break;
  case par_HRUtype::KF_NONLIN:
    value =  pars[26];
    break;
  case par_HRUtype::C:
    value =  pars[27];
    break;
  case par_HRUtype::INFR_MAX:
    value =  pars[28];
    break;
  case par_HRUtype::RF:
    value =  pars[29];
    break;
  case par_HRUtype::WP:
    value =  pars[30];
    break;
  case par_HRUtype::SMAX:
    value =  pars[31];
    break;
  case par_HRUtype::RBAI:
    value =  pars[32];
    break;
  case par_HRUtype::RBEI:
    value =  pars[33];
    break;
}
  return value;

}

/** \brief returns the number of fast runoff reservoirs
 *
 * \return number of fast reservoirs
 *
 */
unsigned params::g_numFastRes() {

  return numFastRes;

}

/** \brief Setting of number of fast reservoirs for fast runoff response
 *
 * \param number of fast reservoirs for fast response
 *
 */
void params::s_numFastRes(const unsigned& numRes) {

  numFastRes = numRes;

  return ;

}

/** \brief Loading the values of to pars sparam type
 *
 * \param vector of pair, each param has value numberSel and ID using par_HRUtype type
 *
 */
void params::s_parLoadToCalib(std::vector<std::pair<numberSel,par_HRUtype>>& parsToLoad) {

  for (auto &  itr : parsToLoad) {
    s_params(itr);
  }

}

/** \brief Loading the default values of HRUunit model
 *
 *
 */
void params::s_default() {

  pars[0] = 2.0;//!< Parameter controlling shape of Pareto distribution of soil storages [0,inf] however [0.5,3],VC1
  pars[1] = 10.0;//!< Max storage of storages distributed by Pareto distribution [0,inf],VC1
  pars[2] = 1.0;//!< Parameter controlling soil evapotranspiration [0,infty] how ever [0.5,3],VC1
  numberSel helpSmax;
  helpSmax = pars[1] / (pars[0]+1);
  pars[3] = helpSmax;//!< Smaxpdm calculate using Cmax and b_soil
  pars[4] = 0.1;//!< Storage coefficient of groundwater storage [0,1],VC1
  pars[5] = 0.5;//!< Storage coefficient of runoff response reservoirs [0,1],VC1
  pars[6] = 0.5;//!< Divider of percolation into the direct flow input par[6]*Perc and  groundwater input (1-par[6])*Perc [0,1],VC1
  pars[7] = 0.02;//!< Divider of gross rainfall as a Canopy input [0,1],VC1
  pars[8] = 0.02;//!< Divider of gross rainfall as a Trunk input [0,1],VC1
  pars[9] = 1.0;//!< The Max canopy storage [0,inf],VC1
  pars[10] = 1.0;//!< The Max stem and trunk storage [0,inf],VC1
  pars[11] = 0.5;//!< The divider of canopy outflow to throughflow and stemflow storage [0,1],VC1
  pars[12] = 1;//!< The threshold temperature for determining snow [-inf,inf] better [-5,5]
  pars[13] = 2;//!< The day degree model for snow melt [o, inf] better [0,2],VC1
  pars[14] = 0;//!<  The threshold temperature for determining melting process [-inf, inf] better [-5,5]
  pars[15] = 1;//!< The maximum capacity of surface retention [0, inf],VC1
  pars[16] = 1.0;//!< Leakage coefficient for linear reservoirs
  pars[17] = 1.0;//!< Direct by-pass coefficient for linear reservoirs
  pars[18] = 1;//!< Power coefficient
  pars[19] = 0.1;//!< Second storage coefficient of groundwater storage [0,1],VC1
  pars[20] = 1.0;//!< THR Threshold coefficient for threshold-controlled linear storage [1,inf]
  pars[21] = 0.5;//!< Divider for two parallel linear reservoirs
  pars[22] = 0;//!< The Cmin for pdm reservoir
  pars[23] = 10;//!< Field capacity
  pars[24] = 0.5;//!< Forest fraction
  pars[25] = 0.5;//!< KF2 Storage coefficient of runoff response reservoirs [0,1]
  pars[26] = 0.5;//!< KF_NONLIN non-linearity parameter, [-] [0,inf)
  pars[27] = 10;//!< C capillary rise [mm d-1] [0, inf)
  pars[28] = 10;//!< INFR_MAX Maximum infiltration rate [mm/d], [0,inf)
  pars[29] = 0.5;//!< RF evaporation reduction factor [-] [0,1]
  pars[30] = 0.5;//!< WP wilting point [-] [0,1]
  pars[31] = 10.0;//!< SMAX soil storage capacity [mm] [0,inf]
  pars[32] = 0.05;//!< RBAI River bank infiltration rate (infiltration to Soil storage)
  pars[33] = 0.1;//!< RBEI River bed infiltration rate (infiltration to groundwater storage)

  numFastRes = 1;

//  std::cout << "Params are initialized." << std::endl;
}

/** \brief Printing the param values to cout
 *
 *
 */
void params::p_param() {

  std::vector<std::string> parsNames {"B_SOIL: ", "C_MAX: ", "B_EVAP: ", "SMAXpdm: ", "KS: ", "KF: ",             \
                                      "ADIV: ", "CDIV: ", "SDIV: ", "CAN_ST: ", "STEM_ST: ", "CSDIV: ", "TETR: ", \
                                      "DDFA: ", "TMEL: ", "RETCAP: ", "L: ", "D_BYPASS: ", "B_EXP: ", "KS2: ",    \
                                      "THR: ", "ALPHA: ","CMIN: ","FC: ","FOREST_FRACT: ", "KF2: ",               \
                                      "KF_NONLIN: ", "C: ", "INFR_MAX: ", "RF: ", "WP: ", "SMAX: ", "RBAI:", "RBEI:"};

  std::cout << std::endl << "Printing the values of parameters:" << std::endl << std::endl;
  for(unsigned pp=0; pp<numPars ; pp++ ) {
    std::cout << parsNames[pp] << pars[pp] << std::endl;
  }
  std::cout << "Number of fast reservoirs "  << numFastRes;
  std::cout << std::endl;
//
  return ;

}

/** \brief Getting the number of params
 *
 *
 */
unsigned params::g_numPars() {

  return(numPars);

}

numberSel params::g_par_low(const par_HRUtype& _parType) {

  numberSel value = 0.0;

  switch(_parType) {
  case par_HRUtype::B_SOIL:
    value = low_pars[0];
    break;
  case par_HRUtype::C_MAX:
    value =  low_pars[1];
    break;
  case par_HRUtype::B_EVAP:
    value =  low_pars[2];
    break;
  case par_HRUtype::SMAXpdm:
    value =  low_pars[3];
    break;
  case par_HRUtype::KS:
    value =  low_pars[4];
    break;
  case par_HRUtype::KF:
    value =  low_pars[5];
    break;
  case par_HRUtype::ADIV:
    value =  low_pars[6];
    break;
  case par_HRUtype::CDIV:
    value =  low_pars[7];
    break;
  case par_HRUtype::SDIV:
    value =  low_pars[8];
    break;
  case par_HRUtype::CAN_ST:
    value =  low_pars[9];
    break;
  case par_HRUtype::STEM_ST:
    value =  low_pars[10];
    break;
  case par_HRUtype::CSDIV:
    value =  low_pars[11];
    break;
  case par_HRUtype::TETR:
    value =  low_pars[12];
    break;
  case par_HRUtype::DDFA:
    value =  low_pars[13];
    break;
  case par_HRUtype::TMEL:
    value =  low_pars[14];
    break;
  case par_HRUtype::RETCAP:
    value =  low_pars[15];
    break;
  case par_HRUtype::L:
    value =  low_pars[16];
    break;
  case par_HRUtype::D_BYPASS:
    value =  low_pars[17];
    break;
  case par_HRUtype::B_EXP:
    value =  low_pars[18];
    break;
  case par_HRUtype::KS2:
    value =  low_pars[19];
    break;
  case par_HRUtype::THR:
    value =  low_pars[20];
    break;
  case par_HRUtype::ALPHA:
    value =  low_pars[21];
    break;
  case par_HRUtype::CMIN:
    value =  low_pars[22];
    break;
  case par_HRUtype::FC:
    value =  low_pars[23];
    break;
  case par_HRUtype::FOREST_FRACT:
    value =  low_pars[24];
    break;
  case par_HRUtype::KF2:
    value =  low_pars[25];
    break;
  case par_HRUtype::KF_NONLIN:
    value =  low_pars[26];
    break;
  case par_HRUtype::C:
    value =  low_pars[27];
    break;
  case par_HRUtype::INFR_MAX:
    value =  low_pars[28];
    break;
  case par_HRUtype::RF:
    value =  low_pars[29];
    break;
  case par_HRUtype::WP:
    value =  low_pars[30];
    break;
  case par_HRUtype::SMAX:
    value =  low_pars[31];
    break;
  case par_HRUtype::RBAI:
    value =  low_pars[32];
    break;
  case par_HRUtype::RBEI:
    value =  low_pars[33];
    break;
  }
  return value;

}

numberSel params::g_par_up(const par_HRUtype& _parType) {

  numberSel value = 0.0;

  switch(_parType) {
  case par_HRUtype::B_SOIL:
    value = up_pars[0];
    break;
  case par_HRUtype::C_MAX:
    value =  up_pars[1];
    break;
  case par_HRUtype::B_EVAP:
    value =  up_pars[2];
    break;
  case par_HRUtype::SMAXpdm:
    value =  up_pars[3];
    break;
  case par_HRUtype::KS:
    value =  up_pars[4];
    break;
  case par_HRUtype::KF:
    value =  up_pars[5];
    break;
  case par_HRUtype::ADIV:
    value =  up_pars[6];
    break;
  case par_HRUtype::CDIV:
    value =  up_pars[7];
    break;
  case par_HRUtype::SDIV:
    value =  up_pars[8];
    break;
  case par_HRUtype::CAN_ST:
    value =  up_pars[9];
    break;
  case par_HRUtype::STEM_ST:
    value =  up_pars[10];
    break;
  case par_HRUtype::CSDIV:
    value =  up_pars[11];
    break;
  case par_HRUtype::TETR:
    value =  up_pars[12];
    break;
  case par_HRUtype::DDFA:
    value =  up_pars[13];
    break;
  case par_HRUtype::TMEL:
    value =  up_pars[14];
    break;
  case par_HRUtype::RETCAP:
    value =  up_pars[15];
    break;
  case par_HRUtype::L:
    value =  up_pars[16];
    break;
  case par_HRUtype::D_BYPASS:
    value =  up_pars[17];
    break;
  case par_HRUtype::B_EXP:
    value =  up_pars[18];
    break;
  case par_HRUtype::KS2:
    value =  up_pars[19];
    break;
  case par_HRUtype::THR:
    value =  up_pars[20];
    break;
  case par_HRUtype::ALPHA:
    value =  up_pars[21];
    break;
  case par_HRUtype::CMIN:
    value =  up_pars[22];
    break;
  case par_HRUtype::FC:
    value =  up_pars[23];
    break;
  case par_HRUtype::FOREST_FRACT:
    value =  up_pars[24];
    break;
  case par_HRUtype::KF2:
    value =  up_pars[25];
    break;
  case par_HRUtype::KF_NONLIN:
    value =  up_pars[26];
    break;
  case par_HRUtype::C:
    value =  up_pars[27];
    break;
  case par_HRUtype::INFR_MAX:
    value =  up_pars[28];
    break;
  case par_HRUtype::RF:
    value =  up_pars[29];
    break;
  case par_HRUtype::WP:
    value =  up_pars[30];
    break;
  case par_HRUtype::SMAX:
    value =  up_pars[31];
    break;
  case par_HRUtype::RBAI:
    value =  up_pars[32];
    break;
  case par_HRUtype::RBEI:
    value =  up_pars[33];
    break;
  }
  return value;

}
