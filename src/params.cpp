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
  numPars = 16;//!< Since the Smax is defined by the Cmax and Bsoil the real number of parameters is numPars-1

  pars.resize(numPars,numPars);
  up_pars.resize(numPars,numPars);
  low_pars.resize(numPars,numPars);
  pars[0] = 2.0;//!< B_SOIL Parameter controlling shape of Pareto distribution of soil storages [0,inf] however [0.5,3],VC1
  pars[1] = 10.0;//!< C_MAX Max storage of storages distributed by Pareto distribution [0,inf],VC1
  pars[2] = 1.0;//!< B_EVAP Parameter controlling soil evapotranspiration [0,infty] how ever [0.5,3],VC1
  numberSel helpSmax;
  helpSmax = pars[1] / (pars[0]+1);
  pars[3] = helpSmax;//!< SMAX Max soil storage calculate using Cmax and b_soil
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
// Upper bounds of parameters
  up_pars[0] = 3.0;//!< B_SOIL Parameter controlling shape of Pareto distribution of soil storages [0,inf] however [0.5,3],VC1
  up_pars[1] = 500.0;//!< C_MAX Max storage of storages distributed by Pareto distribution [0,inf],VC1
  up_pars[2] = 3.0;//!< B_EVAP Parameter controlling soil evapotranspiration [0,infty] how ever [0.5,3],VC1
  numberSel helpSmaxUp;
  helpSmaxUp = up_pars[1] / (up_pars[0]+1);
  up_pars[3] = helpSmaxUp;//!< SMAX Max soil storage calculate using Cmax and b_soil
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
// Lower bounds of parameters
  low_pars[0] = 0.0;//!< B_SOIL Parameter controlling shape of Pareto distribution of soil storages [0,inf] however [0.5,3],VC1
  low_pars[1] = 0.0;//!< C_MAX Max storage of storages distributed by Pareto distribution [0,inf],VC1
  low_pars[2] = 0.0;//!< B_EVAP Parameter controlling soil evapotranspiration [0,infty] how ever [0.5,3],VC1
  numberSel helpSmaxLow;
  helpSmaxLow = low_pars[1] / (low_pars[0]+1);
  low_pars[3] = helpSmaxLow;//!< SMAX Max soil storage calculate using Cmax and b_soil
  low_pars[4] = 0.0;//!< KS Storage coefficient of groundwater storage [0,1],VC1
  low_pars[5] = 0.0;//!< KF Storage coefficient of runoff response reservoirs [0,1],VC1
  low_pars[6] = 0.0;//!< ADIV Divider of percolation into the overflow input par[6]*Perc and  groundwater input (1-par[6])*Perc [0,1],VC1
  low_pars[7] = 0.0;//!< CDIV Divider of gross rainfall as a Canopy input [0,1],VC1
  low_pars[8] = 0.0;//!< SDIV Divider of gross rainfall as a Trunk input [0,1],VC1
  low_pars[9] = 0.0;//!< CAN_ST The Max canopy storage [0,inf],VC1
  low_pars[10] = 0.0;//!< STEM_ST The Max stem and trunk storage [0,inf],VC1
  low_pars[11] = 0.0;//!< CSDIV The divider of canopy outflow to throughflow and stemflow storage [0,1],VC1
  low_pars[12] = 1.0;//!< TETR The threshold temperature for determining snow [-inf,inf] better [-5,5]
  low_pars[13] = 0.0;//!<  DDFA The day degree model for snow melt [o, inf] better [0,2],VC1
  low_pars[14] = 0.0;//!<  TMEL The threshold temperature for determining melting process [-inf, inf] better [-5,5]
  low_pars[15] = 0.0;//!< RETCAP The maximum capacity of surface retention [0, inf],VC1
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
    pars[3] = pars[1] / (pars[0]+1);
//    std::cout << "New b_soil --> loaded\n";
    break;
  case par_HRUtype::C_MAX:
    pars[1] = par_dta;
    pars[3] = pars[1] / (pars[0]+1);
//    std::cout << "New c_max --> loaded\n";
    break;
  case par_HRUtype::B_EVAP:
    pars[2] = par_dta;
//    std::cout << "New b_evap --> loaded\n";
    break;
  case par_HRUtype::SMAX:
    pars[3] = par_dta;
    pars[3] = pars[1] / (pars[0]+1);//For preventing the losing the link between b_soil and cmax
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
  }
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
    pars[3] = pars[1] / (pars[0]+1);
//    std::cout << "New b_soil --> loaded\n";
    break;
  case par_HRUtype::C_MAX:
    pars[1] = par_dta;
    pars[3] = pars[1] / (pars[0]+1);
//    std::cout << "New c_max --> loaded\n";
    break;
  case par_HRUtype::B_EVAP:
    pars[2] = par_dta;
//    std::cout << "New b_evap --> loaded\n";
    break;
  case par_HRUtype::SMAX:
    pars[3] = par_dta;
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
  case par_HRUtype::SMAX:
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
  pars[3] = helpSmax;//!< Smax calculate using Cmax and b_soil
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

  numFastRes = 1;

//  std::cout << "Params are initialized." << std::endl;
}

/** \brief Printing the param values to cout
 *
 *
 */
void params::p_param() {

  std::vector<std::string> parsNames {"B_SOIL: ", "C_MAX: ", "B_EVAP: ", "SMAX: ", "KS: ", "KF: ", \
                                      "ADIV: ", "CDIV: ", "SDIV: ", "CAN_ST: ", "STEM_ST: ", "CSDIV: ", "TETR: ", "DDFA: ", "TMEL: ", "RETCAP: "};

  std::cout << std::endl << "Printing the values of parameters:" << std::endl << std::endl;
  for(unsigned pp=0; pp<numPars ; pp++ ) {
    std::cout << parsNames[pp] << pars[pp] << std::endl;
  }
  std::cout << "Number of fast reservoirs "  << numFastRes;
  std::cout << std::endl;
//
  return ;

}

