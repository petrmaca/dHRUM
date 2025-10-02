#include "dam.h"

/** \brief constructor of dam
 *
 * initialization of dam
 *
 */
dam::dam(): tstRM(0),
hyd_dta(),
prev_DamS(0.0),
damMax(0.0),
MRF(0.0),
damArea(0.0),
damLeng(0.0),
damBank(0.0),
damEt_TYPE{},
damSperc_TYPE{},
damGperc_TYPE{},
damRout_TYPE{}
{
  prev_DamS = get_initState(init_dStype::DAMS);
  tstRM = 0;
}

/** \brief destructor of dam
 *
 *  dam destructor
 *
 */
dam::~dam() {
  //dtor
}


/** \brief copy constructor of dam
 *
 * \param
 * \param
 * \return
 *
 */
dam::dam(const dam& other): tstRM(0),
hyd_dta(),
prev_DamS(0.0),
damMax(0.0),
MRF(0.0),
damArea(0.0),
damLeng(0.0),
damBank(0.0),
damEt_TYPE{},
damSperc_TYPE{},
damGperc_TYPE{},
damRout_TYPE{}
{
  tstRM = other.tstRM;//!< The counter for main loop in run model
  hyd_dta = other.hyd_dta;
  prev_DamS= other.prev_DamS;
  damMax= other.damMax;
  MRF= other.MRF;
  damArea= other.damArea;
  damLeng= other.damLeng;
  damBank= other.damBank;
  damEt_TYPE= other.damEt_TYPE;
  damSperc_TYPE= other.damSperc_TYPE;
  damGperc_TYPE= other.damGperc_TYPE;
  damRout_TYPE= other.damRout_TYPE;
}


/** \brief assignment operator for single pdm unit
 *
 * \param single pdm unit of DHM to be assigned
 *
 * \return single pdm unit
 *
 */
dam& dam::operator=(const dam& rhs) {


  if (this == &rhs) return *this;
  else {
    tstRM = rhs.tstRM;//!< The counter for main loop in run model
    hyd_dta = rhs.hyd_dta;
    prev_DamS= rhs.prev_DamS;
    damMax= rhs.damMax;
    MRF= rhs.MRF;
    damArea= rhs.damArea;
    damLeng= rhs.damLeng;
    damBank= rhs.damBank;
    damEt_TYPE= rhs.damEt_TYPE;
    damSperc_TYPE= rhs.damSperc_TYPE;
    damGperc_TYPE= rhs.damGperc_TYPE;
    damRout_TYPE= rhs.damRout_TYPE;

  } // handle self assignment
  //assignment operator
  return *this;
}



/** \brief initialization of all hdata of all state variables and fluxes using the value val
 *
 * \param value of initialization
 *
 */
void dam::init_inputs(numberSel val, unsigned numDTA) {

  hdata dta(val, numDTA);
  // std::cout << "Initializing the hdata setup " << get_numdta() << " ups numDTA "<< numDTA << std::endl;
  //input time series
  set_data(dta,dam_ts::PREC);
  set_data(dta,dam_ts::DAMS);
  set_data(dta,dam_ts::INFL);
  set_data(dta,dam_ts::DAIS);
  set_data(dta,dam_ts::DAIG);
  set_data(dta,dam_ts::OUFL);
  set_data(dta,dam_ts::ETDM);
  set_data(dta,dam_ts::OFLW);
  set_data(dta,dam_ts::OULT);
  set_data(dta,dam_ts::INLT);
  // std::cout << " tor ok\n";

  // get_numdta();

  // std::cout << "The total number of initialized time intervals is " << get_numdta() << " ." << std::endl;

  return ;

}

/** \brief Sets the inital values of states for all reservoirs to zero
 *
 * \param number of reservoirs in fast response
 *
 */
void dam::set_ZeroinitStates(const unsigned& numres) {
  hdata help_data(0.0,numres);
  numberSel zeroState = 0.0;
  //  help_data;
  hyd_dta.s_initStates(help_data,zeroState,init_dStype::DAMS);
  return ;
}

numberSel dam::get_initState(const init_dStype& _Stype) {

  return hyd_dta.g_initState(_Stype);

}

/** \brief loading data into data vectors
 *
 */
void dam::set_data(const hdata& dta,const dam_ts& _tsType) {

  hyd_dta.s_data(dta,_tsType,true);

  //  std::cout << "Well done data loaded.\n " <<std::endl;
  return ;
}
