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
  set_data(dta,dam_ts::TEMP);
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

numberSel dam::dam_ET(ETdam_type _etdam_type) {
  numberSel Etdam = 0.0;

  switch(_etdam_type) {
  case ETdam_type::ETdam1: {
    //BERAN, A., KAŠPÁREK, L., VIZINA, A. a ŠUHÁJKOVÁ, P. Ztráta vody výparem z volné vodní hladiny. Vodohospodářské technicko-ekonomické informace, 2019, roč. 61, č. 4, str. 12–18. ISSN 0322-8916.
    Etdam = 0.0824 * std::pow(get_dta(tstRM, dam_ts::TEMP),1.289);
    //std::cout<<"ETpond1"<<std::endl;
    break;
  }

  case ETdam_type::ETdam2: {
    //Water surface evapotranspiration [mm/den] (ČSN 75241(str41) prumer je cca 750 mm/rok) ale meni se to hodně v zavislosti na rocnim obdobi v lete až 5 mm/den
    Etdam = 2; //here is only a constant, need another equation
    //std::cout<<"ETpond2"<<std::endl;
    break;
  }
  }

  return Etdam;
}

numberSel dam::Dam_SOISperc(DamSOISPerc_type _soisDam_type) {
  numberSel PoiS = 0.0;

  switch(_soisDam_type) {
  case DamSOISPerc_type::noDamSOISPerc: {
    PoiS = 0.0;
    //std::cout<<"noSOISPerc"<<std::endl;
    break;
  }
  case DamSOISPerc_type::DamSOISPerc1: {
    PoiS = 0.0001; // k [m/s] - coarse sand / gravel
    //std::cout<<"SOISPerc - coarse sand / gravel"<<std::endl;
    break;
  }
  case DamSOISPerc_type::DamSOISPerc2: {
    PoiS = 0.000001; // k [m/s] - sandy loam
    //std::cout<<"SOISPerc - sandy loam"<<std::endl;
    break;
  }
  case DamSOISPerc_type::DamSOISPerc3: {
    PoiS = 0.000000001; // k [m/s] - clay
    //std::cout<<"SOISPerc - clay"<<std::endl;
    break;
  }
  }

  return PoiS;
}

numberSel dam::Dam_GWperc(DamGWPerc_type _gwDam_type) {
  numberSel PoiG = 0.0;

  switch(_gwDam_type) {
  case DamGWPerc_type::noDamGWPerc: {
    PoiG =  0.0;
    //std::cout<<"noGWPerc"<<std::endl;
    break;
  }
  case DamGWPerc_type::DamGWPerc1: {
    PoiG =  0.0001; // k [m/s] - coarse sand / gravel
    //std::cout<<"GWPerc - coarse sand / gravel"<<std::endl;
    break;
  }
  case DamGWPerc_type::DamGWPerc2: {
    PoiG = 0.000001; // k [m/s] - sandy loam
    //std::cout<<"GWPerc - sandy loam"<<std::endl;
    break;
  }
  case DamGWPerc_type::DamGWPerc3: {
    PoiG = 0.000000001; // k [m/s] - clay
    //std::cout<<"GWPerc - clay"<<std::endl;
    break;
  }
  }

  return PoiG;
}

numberSel dam::pond_regular_out(DamRouT_type _RouT_type) {
  numberSel RouT = 0.0; //pond regulat outflow

  switch(_RouT_type) {
  case DamRouT_type::noDamRouT: {
    RouT = 0.0;
    //std::cout<<"noPondRouT"<<std::endl;
    break;
  }
  case DamRouT_type::DamRouT1: {
    //monk with boards - Basin
    numberSel h = 0.1; // [m]
    numberSel b = 0.8; // [m]
    numberSel m = 0.385; //
    RouT= m*b*pow(2*9.81,0.5)*pow(h,3/2);
    //std::cout<<"PondRouT1"<<std::endl;
    break;
  }
  case DamRouT_type::DamRouT2: {
    //pipe - Free outlet through a small hole in the wall
    numberSel mi = 0.8; //hodnota soucinitele
    numberSel d = 0.125; // [m] prumer trubky
    numberSel h = 1; // [m] meni se
    RouT = mi*(3.14*(d*d)/4)*pow((2*h),0.5);
    //std::cout<<"PondRouT2"<<std::endl;
    break;
  }
  case DamRouT_type::DamRouT3: {
    //constant
    RouT = 0.0005; //[m3/s]
    //std::cout<<"PondRouT3"<<std::endl;
    break;
  }
  }

  return RouT;
}

/** \brief Get value of state or flux
 *
 * \param ID of variable
 *
 * \return value of flux or state
 *
 */
numberSel dam::get_dta(const unsigned& tst, const dam_ts& _tsType) {
  return hyd_dta.g_dta(tst, _tsType);
}
