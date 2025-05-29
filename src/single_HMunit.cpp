#include "single_HMunit.h"

/** \brief constructor of single pdm unit
 * initialization of single unit in DHM
 *
 */
single_HMunit::single_HMunit(): tstRM(0),
  par_HRU(),
  hyd_dta(),
  prev_Soil(0.0),
  prev_Grou(0.0),
  prevCanS(0.0),
  prevSteS(0.0),
  prevSnoS(0.0),
  prev_SurS(0.0),
  prev_GroS1(0.0),
  prev_GroS2(0.0),
  et_demand(0.0),
  help_nmbFR(0),
  ifrb(0),
  Area(0),
  IdHru(),
  gs_STORAGE{},
  soil_STORAGE{},
  intrc_STORAGE{}{

  set_nmbFastres(1);
  help_nmbFR = get_nmbFastRes();
//  std::cout << "\nFast runoff response has " << help_nmbFR << " reservoirs." << std::endl;
  set_ZeroinitStates(help_nmbFR);
//  std::cout<< "Soil init " << get_initState(init_Stype::SOIL) << std::endl;
  prev_Soil = get_initState(init_Stype::SOIL);
  prev_Grou = get_initState(init_Stype::GROUNDWAT);
  prevCanS = get_initState(init_Stype::CANS);
  prevSteS  = get_initState(init_Stype::STES);
  prevSnoS = get_initState(init_Stype::SNOS);
  prev_SurS = get_initState(init_Stype::SURFRET);
  prev_GroS1 = get_initState(init_Stype::GROS1);
  prev_GroS2 = get_initState(init_Stype::GROS2);
  p_defaultParams(false);
//  std::cout << "INITprevDR " << prev_Grou << std::endl;
  tstRM = 0;
  gs_STORAGE = gs_STORtype::LIN_RES;
  soil_STORAGE = soil_STORtype::PDM;
  intrc_STORAGE = interception_STORtype::Rutter_Gash;
  et_demand = 0.0;

}

/** \brief destructor of single pdm unit
 *
 *  destructor single unit in DHM
 *
 */
single_HMunit::~single_HMunit() {

  //dtor

}


/** \brief copy constructor of single pdm unit
 *
 * \param
 * \param
 * \return
 *
 */
single_HMunit::single_HMunit(const single_HMunit& other): tstRM(0),
par_HRU(),
hyd_dta(),
prev_Soil(0.0),
prev_Grou(0.0),
prevCanS(0.0),
prevSteS(0.0),
prevSnoS(0.0),
prev_SurS(0.0),
prev_GroS1(0.0),
prev_GroS2(0.0),
et_demand(0.0),
help_nmbFR(0),
ifrb(0),
Area(0),
IdHru(),
gs_STORAGE{},
soil_STORAGE{},
intrc_STORAGE{}
{

  tstRM = other.tstRM;//!< The counter for main loop in run model
  par_HRU = other.par_HRU;//!< The parameters in PDM instances
  hyd_dta = other.hyd_dta;//!< The data of all time series of hydrological variables
  prev_Soil = other.prev_Soil;//!< The helper variable for updating soil storage
  prev_Grou = other.prev_Grou;//!< The helper variable for updating groundwater storage
  prevCanS = other.prevCanS;//!<  The helper variable for Canopy interception storage
  prevSteS = other.prevSteS;//!<  The helper variable for Stem interception storage
  prevSnoS = other.prevSnoS;//!<  The helper variable for Snow storage
  prev_SurS = other.prev_SurS;//!< The helper variable for updating surface storage
  prev_GroS1 = other.prev_GroS1;
  prev_GroS2 = other.prev_GroS2;
  et_demand = other.et_demand;
  help_nmbFR = other.help_nmbFR;//!< The helper for number of fast reservoirs
  ifrb = other.ifrb;//!< For loop counter
  Area = other.Area;//!< The area of HM unit in m2
  IdHru = other.IdHru;
  gs_STORAGE = other.gs_STORAGE;
  soil_STORAGE = other.soil_STORAGE;
  intrc_STORAGE = other.intrc_STORAGE;

}


/** \brief assignment operator for single pdm unit
 *
 * \param single pdm unit of DHM to be assigned
 *
 * \return single pdm unit
 *
 */
single_HMunit& single_HMunit::operator=(const single_HMunit& rhs) {


  if (this == &rhs) return *this;
  else {
    tstRM = rhs.tstRM;//!< The counter for main loop in run model
    par_HRU = rhs.par_HRU;//!< The parameters in PDM instances
    hyd_dta = rhs.hyd_dta;//!< The data of all time series of hydrological variables
    prev_Soil = rhs.prev_Soil;//!< The helper variable for updating soil storage
    prev_Grou = rhs.prev_Grou;//!< The helper variable for updating groundwater storage
    prevCanS = rhs.prevCanS;//!<  The helper variable for Canopy interception storage
    prevSteS = rhs.prevSteS;//!<  The helper variable for Stem interception storage
    prevSnoS = rhs.prevSnoS;//!<  The helper variable for Snow storage
    prev_SurS = rhs.prev_SurS;//!< The helper variable for updating surface storage
    prev_GroS1 = rhs.prev_GroS1;
    prev_GroS2 = rhs.prev_GroS2;
    et_demand = rhs.et_demand;
    help_nmbFR = rhs.help_nmbFR;//!< The helper for number of fast reservoirs
    ifrb = rhs.ifrb;//!< For loop counter
    Area = rhs.Area;//!< The area of HM unit in m2
    IdHru = rhs.IdHru;//!< The ID of HRU
    gs_STORAGE = rhs.gs_STORAGE;//!< The type of groundwater storage
    soil_STORAGE = rhs.soil_STORAGE;
    intrc_STORAGE = rhs.intrc_STORAGE;

  } // handle self assignment
  //assignment operator
  return *this;
}


void single_HMunit::set_varValue(const numberSel& dta,const unsigned& tst,const ts_type& _tsType) {

  hyd_dta.s_varVal(dta, tst, _tsType);

  return ;

}



/** \brief loading data into data vectors
 *
 */
void single_HMunit::set_data(const hdata& dta,const ts_type& _tsType) {

  hyd_dta.s_data(dta,_tsType,true);

  //  std::cout << "Well done data loaded.\n " <<std::endl;
  return ;
}

/** \brief loading precipitation and temperature data only into data valarrays
 *
 */
void single_HMunit::set_data_prec_temp(const hdata& _prec_dta,const hdata& _temp_dta) {

  hyd_dta.s_data(_prec_dta,ts_type::PREC,true);
  hyd_dta.s_data(_temp_dta,ts_type::TEMP,false);

  //  std::cout << "Well done precipitation and temperature data loaded.\n" << std::endl;
  return ;
}

unsigned single_HMunit::get_numdta() {

  return hyd_dta.g_numdta();

}

/** \brief Get value of parameter
 *
 * \param ID of parameter
 *
 * \return value of parameter
 *
 */
numberSel single_HMunit::get_par(const par_HRUtype& _parType) {

  return par_HRU.g_par(_parType);

}

/** \brief Get value of state or flux
 *
 * \param ID of variable
 *
 * \return value of flux or state
 *
 */
numberSel single_HMunit::get_dta(const unsigned& tst, const ts_type& _tsType) {



  return hyd_dta.g_dta(tst, _tsType);

}


void single_HMunit::set_nmbFastres(const unsigned& nmbRes) {

  par_HRU.s_numFastRes(nmbRes);
  hyd_dta.s_numFastRes(nmbRes);

  hdata help_data(0.0,nmbRes);
  hyd_dta.s_initStates(help_data,0.0,init_Stype::FASTRES);
  help_nmbFR = nmbRes;

//
//   for(unsigned i=0;i<nmbRes;i++){
//     std::cout << hyd_dta.g_oneFastResOut(i) << " hh \n ";
//   }

  // std::cout << nmbRes << " hh " << hyd_dta.g_numFastRes() << "\n";
  // std::cout << nmbRes << " pp " << par_HRU.g_numFastRes()<< "\n";

  return;
}

unsigned single_HMunit::get_nmbFastRes() {


  return par_HRU.g_numFastRes();

}

/** \brief Get init states
 *
 * \param ID of state variable
 *
 * \return value of init state
 *
 */
numberSel single_HMunit::get_initState(const init_Stype& _Stype) {

  return hyd_dta.g_initState(_Stype);

}

/** \brief Sets the inital values of states for all reservoirs to zero
 *
 * \param number of reservoirs in fast response
 *
 */
void single_HMunit::set_ZeroinitStates(const unsigned& numres) {

  hdata help_data(0.0,numres);
  numberSel zeroState = 0.0;
  //  help_data;
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::SOIL);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::GROUNDWAT);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::SURFRET);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::CANS);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::STES);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::FASTRES);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::SNOS);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::SURFRET);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::GROS1);
  hyd_dta.s_initStates(help_data,zeroState,init_Stype::GROS2);

  return ;
}

/** \brief get value of output of fast response linear reservoir
 *
 * \param ID of linear reservoir
 *
 */
numberSel single_HMunit::get_outFastRes(const unsigned& itFasRes) {

  return hyd_dta.g_oneFastResOut(itFasRes);

}

/** \brief get value of state of fast response linear reservoir
 *
 * \param ID of linear reservoir
 *
 */
numberSel single_HMunit::get_stateFastres(const unsigned& itFasRes) {

  return hyd_dta.g_oneFastResState(itFasRes);

}

/** \brief Setting of output from reservoirs into the input to next reservoir
 *
 * \param output to be input for next reservoir
 * \param ID of output
 *
 */
void single_HMunit::set_outFastRes(const numberSel& helpOut,const unsigned& itFasRes) {

  hyd_dta.s_oneFastreResUout(helpOut,itFasRes);

  return ;
}

/** \brief Setting of state to the system of fast reservoirs
 *
 * \param state of reservoir with ID
 * \param ID of reservoir
 *
 */

void single_HMunit::set_stateFastRes(const numberSel& helpState,const unsigned& itFasRes) {

  hyd_dta.s_oneFastreResUstate(helpState,itFasRes);

  return ;
}

/** \brief Update all surface retention
 *
 */
void single_HMunit::surface_retention() {

  numberSel RetOut = 0.0, EvapSR = 0.0;
  //  RetOut = std::max((static_cast<numberSel>(prev_SurS) - static_cast<numberSel>(get_par(par_HRUtype::RETCAP))),0.0);
  RetOut = std::max((prev_SurS - get_par(par_HRUtype::RETCAP)),0.0);
  // std::cout << RetOut << "  retout " << tstRM <<std::endl;
  prev_SurS = prev_SurS - RetOut;

  // Evaporation according to Beran VTEI
  // if (get_dta(tstRM, ts_type::TEMP) > 0.0) {
  //   EvapSR = 0.0824 * std::pow(get_dta(tstRM, ts_type::TEMP),1.289);
  // } else EvapSR = 0.0;

  // EvapSR = std::max((static_cast<numberSel>(0.0824 * std::pow(get_dta(tstRM, ts_type::TEMP),1.289))),0.0);
  EvapSR = 0.0824 * std::pow(get_dta(tstRM, ts_type::TEMP),1.289);

  if (EvapSR > prev_SurS) {
    EvapSR = prev_SurS;
  }

  if((EvapSR < 0.0)||(std::isnan(EvapSR))) {
    EvapSR = 0.0;
  }

  numberSel help_EvapSR = update_ETDEMAND(EvapSR, false);
  et_demand = update_ETDEMAND(EvapSR, true);
  EvapSR = help_EvapSR;

  // if(EvapSR >= et_demand) {
  //   EvapSR = et_demand;
  //   et_demand = 0.0;
  // }
  // else {
  //   et_demand = et_demand - EvapSR;
  //   }



   // std::cout << EvapSR << "  EvapSR " << tstRM << " " << get_dta(tstRM, ts_type::TEMP) << " u " << (std::pow((-1.5),1.289)) <<std::endl;
  // if(tstRM == 469) std::cout <<  " e " << EvapSR <<std::endl;
  // prev_SurS = std::max((prev_SurS - EvapSR),0.0);
  prev_SurS = prev_SurS - EvapSR;
  // if(tstRM == 469) std::cout << RetOut << " ps " << prev_SurS << " e " << EvapSR <<std::endl;
  set_varValue(prev_SurS, tstRM, ts_type::SURS);

  if(get_dta(tstRM, ts_type::TEMP) < get_par(par_HRUtype::TETR)) {
    prev_SurS = prev_SurS  + get_dta(tstRM, ts_type::TROF) +  \
      (1 - get_par(par_HRUtype::CDIV)  - get_par(par_HRUtype::SDIV)) * (get_dta(tstRM, ts_type::MELT));
  } else {
    prev_SurS = prev_SurS  + get_dta(tstRM, ts_type::TROF) +  \
      (1 - get_par(par_HRUtype::CDIV)  - get_par(par_HRUtype::SDIV)) * (get_dta(tstRM, ts_type::MELT) + (1 - get_par(par_HRUtype::CDIV)  - get_par(par_HRUtype::SDIV)) *get_dta(tstRM, ts_type::PREC));
      }

  // set_varValue(EvapSR, tstRM, ts_type::AET);
  set_varValue(EvapSR, tstRM, ts_type::ETSW);
  set_varValue(RetOut,tstRM,ts_type::PREF);
  // set_varValue(prev_SurS, tstRM, ts_type::SURS);

  return ;
}


/** \brief updates states in soil buffer in single pdm unit
 *
 *  updating of states in soil buffer in single pdm unit
 *  calculating related fluxes
 *
 */
void single_HMunit::soil_buffer(soil_STORtype _soil_STORtype) {

  numberSel c_init = 0.0, overFl1 = 0.0, ppInf = 0.0, overFl2 = 0.0, overFl3 = 0.0, c_prop = 0.0, next_soil = 0.0, overFL = 0.0, evap = 0.0, aet = 0.0, c_contr=0.0, E_b = 0.0, E_v = 0.0, E_s = 0.0,E_p = 0.0, P_s = 0.0, P_n = 0.0;

switch(_soil_STORtype) {

case soil_STORtype::PDM: {

  numberSel diff = 0.0, diff2= 0.0;
  // Estimation of Soil Water Depth C using total soil basin Storage S from previous day see Moore description of PDM HESS 2007,
  //  Wood and bascics from 1992
  // Eric F. Wood, D. P. Lettenmaier, V. G. Zartarian A land-surface hydrology parameterization with subgrid variability for general circulation models
  //   // eq 3a or 18a Vic paper

  // prev_soil is soil storage from previous time step

  if(prev_Soil <= get_par(par_HRUtype::CMIN)){
// the storages in soil filled with c<=c_min
// no pdm Snew = Sold + PrecEf - evap from c with c<=c_min
    prev_Soil = prev_Soil + get_dta(tstRM, ts_type::PREF);
// checking the ooverflow during one interval
    if(prev_Soil > get_par(par_HRUtype::SMAX)){
      overFL = prev_Soil - get_par(par_HRUtype::SMAX);
      prev_Soil = prev_Soil - overFL;
    }

    evap  = get_dta(tstRM, ts_type::PET)* prev_Soil /  get_par(par_HRUtype::C_MAX);
    numberSel help_EvapSR = update_ETDEMAND(evap, false);
    et_demand = update_ETDEMAND(evap, true);
    evap = help_EvapSR;

  } else{
// c_init height in pdm soil storage
     if((get_par(par_HRUtype::SMAX) - prev_Soil)>=0.0){
      c_init = get_par(par_HRUtype::CMIN) + (get_par(par_HRUtype::C_MAX)-get_par(par_HRUtype::CMIN)) * (1 - pow(((get_par(par_HRUtype::SMAX) - prev_Soil) / (get_par(par_HRUtype::SMAX)- get_par(par_HRUtype::CMIN))),(1/(get_par(par_HRUtype::B_SOIL) + 1))));
     } else {
       c_init = get_par(par_HRUtype::C_MAX);
       diff = prev_Soil - get_par(par_HRUtype::SMAX);
     }

    //Overflow if soil tank fully filled
    overFl1 = std::max(static_cast<numberSel>(get_dta(tstRM, ts_type::PREF) + c_init - get_par(par_HRUtype::C_MAX)), static_cast<numberSel>(0.0));
    if(overFl1 > get_dta(tstRM, ts_type::PREF)) {
      ppInf = 0.0;
    } else {
      ppInf = get_dta(tstRM, ts_type::PREF) - overFl1;
    }

    //Newly proposed soil water depth C
    c_prop = std::min(static_cast<numberSel>(ppInf + c_init), static_cast<numberSel>(get_par(par_HRUtype::C_MAX)));
    // newly proposed storage
    next_soil = get_par(par_HRUtype::CMIN) + (get_par(par_HRUtype::SMAX)-get_par(par_HRUtype::CMIN)) * (1 - pow(( (get_par(par_HRUtype::C_MAX) -c_prop )/ (get_par(par_HRUtype::C_MAX)-get_par(par_HRUtype::CMIN))),(get_par(par_HRUtype::B_SOIL) + 1)));
  //Overflow for small C than C_max according to Jherman https://github.com/jdherman/hymod/blob/master/HyMod.cpp
  overFl2 = std::max(static_cast<numberSel>(0.0),static_cast<numberSel>(ppInf - next_soil + prev_Soil + diff));

  if((get_par(par_HRUtype::SMAX) - next_soil)<0){
    prev_Soil = next_soil;
  } else {
    diff2 =  next_soil - get_par(par_HRUtype::SMAX);
    prev_Soil = get_par(par_HRUtype::SMAX);
  }

  if((get_par(par_HRUtype::SMAX) - prev_Soil)>0) {
      evap =  std::min(prev_Soil,(get_dta(tstRM, ts_type::PET)*(1 - pow(((get_par(par_HRUtype::SMAX) - prev_Soil) / get_par(par_HRUtype::SMAX)), (get_par(par_HRUtype::B_EVAP))))));
  } else evap =  std::min(prev_Soil,(get_dta(tstRM, ts_type::PET)));


  numberSel help_evap = update_ETDEMAND(evap, false);
  et_demand = update_ETDEMAND(evap, true);
  evap = help_evap;

// if(tstRM == 35){
     // std::cout <<" before et demand  "<<  evap << " etdem "<< et_demand << " pet evap " << get_dta(tstRM, ts_type::PET) << " help evap " << help_evap << std::endl;
  // }
  //Soil buffer state
  // numberSel help_EvapSR = update_ETDEMAND(evap, false);
  // et_demand = update_ETDEMAND(evap, true);
  // evap = help_EvapSR;

  // std::cout << " evap 1 " << evap <<  " next_soil " << next_soil <<"\n";

  prev_Soil = std::max(static_cast<numberSel>(prev_Soil - evap),static_cast<numberSel>(0.0));

  // std::cout << "ending evap  " << evap <<  " prev_Soil " << prev_Soil << " first part "<< (get_par(par_HRUtype::SMAX) - prev_Soil) / get_par(par_HRUtype::SMAX) <<"\n";
  //Total overflow
  // std::cout << "diff " << diff << " diff 2 " << diff2 << std::endl; // by the definition of pdm diff + diff2 = 0
  overFL = overFl1 + overFl2 ;
  }

  if(prev_Soil - (get_par(par_HRUtype::SMAX) )>0){
    overFL = overFL + (prev_Soil - (get_par(par_HRUtype::SMAX)));
    prev_Soil =  (get_par(par_HRUtype::SMAX));
  }


  // if(next_soil<0.0) next_soil=0.0;
  // if(std::isnan(next_soil)) next_soil=0.0;
  set_varValue(prev_Soil, tstRM, ts_type::SOIS);
  set_varValue(evap,tstRM, ts_type::EVBS);
  set_varValue(overFL, tstRM, ts_type::PERC);

  // aet = get_dta(tstRM,ts_type::EVAC) +  get_dta(tstRM,ts_type::EVAS)  +  evap;
  // set_varValue(aet, tstRM, ts_type::AET);
  // std::cout <<  " prev_soil 1 " << prev_Soil <<"\n";
  // prev_Soil = next_soil;



  // std::cout <<  " prev_soil 2 " << prev_Soil <<"\n";


  // std::cout << " et2 " << et_demand << " evap2 " << evap << "\n";
  break;
}

// Alberto Montanari implementation from hi R package
// https://github.com/albertomontanari/hymodbluecat/blob/main/src/hymodfortran.f95
case soil_STORtype::PDM2: {
  numberSel w1 = 0.0, dummy = 0.0, c1 = 0.0, c2 = 0.0, er1 =0.0, w2 = 0.0, er2 = 0.0,evap = 0.0, overFL =0.0;

  w1 = prev_Soil;

  dummy = (1 - ((get_par(par_HRUtype::B_SOIL)+1)*w1/get_par(par_HRUtype::C_MAX)));
  dummy = std::max(dummy,0.0);

  c1 = get_par(par_HRUtype::C_MAX) * (1- (std::pow(dummy,(1/(get_par(par_HRUtype::B_SOIL)+1)))));

  c2 = std::min((c1 + get_dta(tstRM, ts_type::PREF)), get_par(par_HRUtype::C_MAX));
  //the  lost of water
  c2 = std::max(c2,0.0);

  er1=std::max(get_dta(tstRM, ts_type::PREF) - get_par(par_HRUtype::C_MAX) +c1,0.0);

  dummy = 1- c2 / get_par(par_HRUtype::C_MAX);
  dummy = std::max(dummy,0.0);

  w2 = (get_par(par_HRUtype::C_MAX) / (get_par(par_HRUtype::B_SOIL)+1)) * (1-(std::pow(dummy,(get_par(par_HRUtype::B_SOIL)+1))));

  er2 = std::max(((c2-c1) - (w2-w1)),0.0);
  evap = (1- (((get_par(par_HRUtype::C_MAX) - c2)/(get_par(par_HRUtype::B_SOIL)+1))/(get_par(par_HRUtype::C_MAX)/(get_par(par_HRUtype::B_SOIL)+1)))) * get_dta(tstRM, ts_type::PET);

  w2 = std::max(w2-evap,0.0);
  //the  lost of et

  overFL = er1 + er2;

  prev_Soil = w2;
  set_varValue(w2, tstRM, ts_type::SOIS);
  set_varValue(evap,tstRM, ts_type::EVBS);
  set_varValue(overFL, tstRM, ts_type::PERC);
  }


case soil_STORtype::COLLIE_V2: {
    //actual evaporation is split between bare soil evaporation E_b and transpiration through vegetation E_v
    E_b = prev_Soil / get_par(par_HRUtype::SMAX) * (1 - get_par(par_HRUtype::FOREST_FRACT)) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET));

    if (prev_Soil > get_par(par_HRUtype::FC)) {
      E_v = get_par(par_HRUtype::FOREST_FRACT) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET));
      //overFl2 [mm/d] is subsurface flow regulated by runoff coefficient KF
      overFl2 = get_par(par_HRUtype::KF) * (prev_Soil - get_par(par_HRUtype::FC));
    } else {
      E_v = get_par(par_HRUtype::FOREST_FRACT) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET)) * (prev_Soil/get_par(par_HRUtype::FC));
      overFl2 = 0;
    }

    if (prev_Soil > get_par(par_HRUtype::SMAX)) {
      //overFl1 [mm/d] is saturation excess overland flow
      overFl1 = static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));
    } else {
      overFl1 = 0;
    }

    next_soil = prev_Soil + static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));

    std::vector<numberSel> updated_vals;

    updated_vals = water_balance(next_soil, evap, std::vector<numberSel>{E_v, E_b});
    next_soil = updated_vals[0];
    evap = updated_vals[1];

    updated_vals.clear();

    updated_vals = water_balance(next_soil, overFL, std::vector<numberSel>{overFl1, overFl2});
    next_soil = updated_vals[0];
    overFL = updated_vals[1];

    updated_vals.clear();

    set_varValue(next_soil, tstRM, ts_type::SOIS);
    set_varValue(evap,tstRM, ts_type::EVBS);
    set_varValue(overFL, tstRM, ts_type::PERC);

    prev_Soil = next_soil;

    break;
}


case soil_STORtype::NEW_ZEALAND: {
    if(prev_Soil > get_par(par_HRUtype::FC)) {
      //E_v is evaporation through vegetation
      E_v = get_par(par_HRUtype::FOREST_FRACT) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET));
    } else {
      E_v = get_par(par_HRUtype::FOREST_FRACT) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET)) * (next_soil / get_par(par_HRUtype::FC));
    }

    //E_b is bare soil evaporation
    E_b = prev_Soil / get_par(par_HRUtype::SMAX) * (1 - get_par(par_HRUtype::FOREST_FRACT)) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET));

    if(prev_Soil >= get_par(par_HRUtype::FC)) {
      //overFl2 is subsurface runoff
      overFl2 = pow(get_par(par_HRUtype::KF) * (prev_Soil - get_par(par_HRUtype::FC)),get_par(par_HRUtype::KF_NONLIN));
    } else {
      overFl2 = 0;
    }

    if(prev_Soil >= get_par(par_HRUtype::SMAX)) {
      //saturation excess runoff
      overFl1 = static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));
    } else {
      overFl1 = 0;
    }

    //overFl3 represents baseflow controlled by time scale parameter KF2
    overFl3 = get_par(par_HRUtype::KF2) * prev_Soil;
    //std::cout << overFl3 << " " << prev_Soil << " " << std::endl;

    next_soil = prev_Soil + static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));

    std::vector<numberSel> updated_vals;

    updated_vals = water_balance(next_soil, evap, std::vector<numberSel>{E_v, E_b});
    next_soil = updated_vals[0];
    evap = updated_vals[1];

    updated_vals.clear();

    updated_vals = water_balance(next_soil, overFL, std::vector<numberSel>{overFl1, overFl2, overFl3});
    next_soil = updated_vals[0];
    overFL = updated_vals[1];

    updated_vals.clear();

    set_varValue(next_soil, tstRM, ts_type::SOIS);
    set_varValue(evap,tstRM, ts_type::EVBS);
    set_varValue(overFL, tstRM, ts_type::PERC);

    prev_Soil = next_soil;

    break;
}

case soil_STORtype::GR4J: {
    //P_n is the net precipitation
    if(static_cast<numberSel>(get_dta(tstRM, ts_type::PREF)) >= static_cast<numberSel>(get_dta(tstRM, ts_type::PET))) {
      P_n = static_cast<numberSel>(get_dta(tstRM, ts_type::PREF)) - static_cast<numberSel>(get_dta(tstRM, ts_type::PET));
    } else {
      P_n = 0;
    }

    //P_s is fraction of the net precipitation P_n redirected to soil moisture
    P_s = P_n * (1 - std::pow(prev_Soil / get_par(par_HRUtype::SMAX), 2));
    //P_s = static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));

    //evap is the net evaporation
    if (static_cast<numberSel>(get_dta(tstRM, ts_type::PET)) > static_cast<numberSel>(get_dta(tstRM, ts_type::PREF))) {
      E_p = static_cast<numberSel>(get_dta(tstRM, ts_type::PET)) - static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));
    } else {
      E_p = 0;
    }

    //E_s is fraction of the net evaporation evap subtracted from soil moisture
    E_s = E_p * (2 * (prev_Soil/get_par(par_HRUtype::SMAX)) - std::pow(prev_Soil/get_par(par_HRUtype::SMAX), 2));

    //overFl is the percolation to deeper soil waters
    overFl1 = std::pow(get_par(par_HRUtype::SMAX), -4) / 4 * (std::pow(4.0/9.0, -4)) * std::pow(std::abs(prev_Soil),5);

    next_soil = prev_Soil + P_s;

    std::vector<numberSel> updated_vals;

    updated_vals = water_balance(next_soil, evap, std::vector<numberSel>{E_s});
    next_soil = updated_vals[0];
    evap = updated_vals[1];

    updated_vals.clear();

    updated_vals = water_balance(next_soil, overFL, std::vector<numberSel>{overFl1});
    next_soil = updated_vals[0];
    overFL = updated_vals[1];

    updated_vals.clear();

    set_varValue(next_soil, tstRM, ts_type::SOIS);
    set_varValue(evap,tstRM, ts_type::EVBS);
    set_varValue(overFL, tstRM, ts_type::PERC);

    prev_Soil = next_soil;

    break;
}

case soil_STORtype::SBROOK_V1: {
    //E_b is bare soil evaporation
    E_b = prev_Soil / get_par(par_HRUtype::SMAX) * (1 - get_par(par_HRUtype::FOREST_FRACT)) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET));

    //E_v is transpiration from vegetation
    //overFl2 is non-linear subsurface flow using the wilting point FC as
    //a threshold for flow generation and two flow parameters KF[d] and RUONOFF_NONLIN[-]
    if(prev_Soil > get_par(par_HRUtype::FC)) {
      E_v = get_par(par_HRUtype::FOREST_FRACT) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET));
      overFl2 = std::pow((prev_Soil - get_par(par_HRUtype::FC)) / get_par(par_HRUtype::KF), 1/get_par(par_HRUtype::KF_NONLIN));
    } else {
      E_v = prev_Soil / get_par(par_HRUtype::FC) * get_par(par_HRUtype::FOREST_FRACT) * static_cast<numberSel>(get_dta(tstRM, ts_type::PET));
      overFl2 = 0;
    }

    //overFl1 is saturation excess from flow [mm/d]

    if(prev_Soil >= get_par(par_HRUtype::SMAX)) {
      overFl1 = static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));
    } else {
      overFl1 = 0;
    }

    next_soil = prev_Soil + static_cast<numberSel>(get_dta(tstRM, ts_type::PREF));

    std::vector<numberSel> updated_vals;

    updated_vals = water_balance(next_soil, evap, std::vector<numberSel>{E_v, E_b});
    next_soil = updated_vals[0];
    evap = updated_vals[1];

    updated_vals.clear();

    updated_vals = water_balance(next_soil, overFL, std::vector<numberSel>{overFl1, overFl2});
    next_soil = updated_vals[0];
    overFL = updated_vals[1];

    updated_vals.clear();

    set_varValue(next_soil, tstRM, ts_type::SOIS);
    set_varValue(evap,tstRM, ts_type::EVBS);
    set_varValue(overFL, tstRM, ts_type::PERC);

    prev_Soil = next_soil;

    break;
}

case soil_STORtype::HILLSLOPE: {

  numberSel perc = 0.0;

    // //OLD 1st version
    // //incomming interception is reduced by interception storage which is assumed to evaporate before the next
    // //precipitation event
    // //P_n = std::max(static_cast<numberSel>(get_dta(tstRM, ts_type::PREF)) - static_cast<numberSel>(get_dta(tstRM, ts_type::INTS)),static_cast<numberSel>(0.0));
    //
    // //Evaporation from soil moisture evap [mm/d] occurs at the potential rate PET whenever possible
    // if(prev_Soil > 0) {
    //   E_b = static_cast<numberSel>(get_dta(tstRM, ts_type::PET));
    // } else {
    //   E_b = 0;
    // }
    //
    // //Storage excess surface runoff overFL [mm/d] depends on the fraction of the catchment that is currently saturated,
    // //expressed through parameters SMAX [mm] and KF_NONLIN
    // // overFl1 = (1 - std::pow(1 - prev_Soil / get_par(par_HRUtype::SMAX), get_par(par_HRUtype::KF_NONLIN))) * P_n;
    //
    //
    // next_soil = prev_Soil + P_n + get_par(par_HRUtype::C) * overFl1;
    //
    // overFl1 = overFl1 * (1-get_par(par_HRUtype::C));
    //
    // std::vector<numberSel> updated_vals;
    //
    // updated_vals = water_balance(next_soil, evap, std::vector<numberSel>{E_b});
    // next_soil = updated_vals[0];
    // evap = updated_vals[1];
    //
    // updated_vals.clear();
    //
    // updated_vals = water_balance(next_soil, overFL, std::vector<numberSel>{overFl1});
    // next_soil = updated_vals[0];
    // overFL = updated_vals[1];
    //
    // updated_vals.clear();
    //
    // set_varValue(next_soil, tstRM, ts_type::SOIS);
    // set_varValue(evap,tstRM, ts_type::EVBS);
    // set_varValue(overFL, tstRM, ts_type::PERC);
    //
    // prev_Soil = next_soil;

    //Evaporation from soil moisture evap [mm/d] occurs at the potential rate PET whenever possible
    //PET calcculated without consideration on et demand
    if(prev_Soil > 0) {
      E_b = static_cast<numberSel>(get_dta(tstRM, ts_type::PET));
    } else {
      E_b = 0;
    }

    if(prev_Soil > E_b) {
      prev_Soil = prev_Soil - E_b;
    } else {
      E_b = prev_Soil;
      prev_Soil = 0.0;
    }

    // Percolation and soil storage update
    // perc = (1-(pow((1 - prev_Soil / get_par(par_HRUtype::SMAX)),get_par(par_HRUtype::KF_NONLIN)))) * get_dta(tstRM,ts_type::PREF);

    if(((prev_Soil + get_dta(tstRM,ts_type::PREF))) > get_par(par_HRUtype::SMAX)) {
     // prev_Soil = get_par(par_HRUtype::SMAX);
      perc = prev_Soil + get_dta(tstRM,ts_type::PREF) - get_par(par_HRUtype::SMAX);
    } else {
      perc = (1-(pow((1 - prev_Soil / get_par(par_HRUtype::SMAX)),get_par(par_HRUtype::KF_NONLIN)))) * get_dta(tstRM,ts_type::PREF);
      }

    prev_Soil = prev_Soil + get_dta(tstRM,ts_type::PREF) - perc;

    set_varValue(prev_Soil, tstRM, ts_type::SOIS);
    set_varValue(E_b,tstRM, ts_type::EVBS);
    set_varValue(perc, tstRM, ts_type::PERC);

    break;
}

case soil_STORtype::PLATEAU: {
    //ppInf is infiltration based on maximum infiltration rate INFR_MAX [mm/d] and infiltration excess [mm/d]
    ppInf = std::min(static_cast<numberSel>(get_dta(tstRM, ts_type::PREF)), get_par(par_HRUtype::INFR_MAX));

    //c_init represents the wilting point
    c_init = get_par(par_HRUtype::WP) * get_par(par_HRUtype::SMAX);

    //evac is evaporation from soil moisture [mm/d] which occurs at the potential rate PET when
    //prev_Soil is above wilting point c_init [mm] and is further constrained by coefficient RF [-]
    evap = static_cast<numberSel>(get_dta(tstRM, ts_type::PET)) * std::max((get_par(par_HRUtype::RF) * (prev_Soil - c_init)/(get_par(par_HRUtype::SMAX) - c_init)), static_cast<numberSel>(0.0));

    //overFL is storage excess
    if(prev_Soil == get_par(par_HRUtype::SMAX)) {
      overFL = ppInf + get_par(par_HRUtype::C);
    } else {
      overFL = 0;
    }

    next_soil = prev_Soil + ppInf + get_par(par_HRUtype::C) * overFL;

    overFL = overFL * (1-get_par(par_HRUtype::C) );

    if(next_soil >= 0) {
      if(next_soil < evap) {
        evap = next_soil;
      }
      next_soil = next_soil - evap;

      if(next_soil >= 0) {
        if(next_soil < overFL) {
          overFL = next_soil;
        }
        next_soil = next_soil - overFL;
      }

    }

    set_varValue(next_soil, tstRM, ts_type::SOIS);
    set_varValue(evap,tstRM, ts_type::EVBS);
    set_varValue(overFL, tstRM, ts_type::PERC);

    prev_Soil = next_soil;

    break;
  }
}
  return ;
}

/** \brief updates states in groundwater storage in single pdm unit
 *
 *  updating of states in groundwater storage in the single pdm unit
 *  calculating related fluxes, groundwater storage represented by a linear reservoir
 *
 */
/*void single_HMunit::slow_response() {

  numberSel BaseOut = 0.0;
  // numberSel Evap = 0.0;



  BaseOut = prev_Grou * get_par(par_HRUtype::KS) ;
  // Evap = prev_Grou * get_par(par_HRUtype::KS) ;
  prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;//- Evap

  set_varValue(BaseOut, tstRM, ts_type::BASF);
  set_varValue(prev_Grou, tstRM,ts_type::GROS);

  return ;
}*/

void single_HMunit::slow_response(gs_STORtype _gs_STORtype) {

  numberSel BaseOut = 0.0;
  numberSel BaseOut_1 = 0.0;
  numberSel BaseOut_2 = 0.0;
  // numberSel prev_Grou1 = 0.0;
  // numberSel prev_Grou2 = 0.0;

  switch(_gs_STORtype) {

  case gs_STORtype::LIN_RES:
    BaseOut = prev_Grou * get_par(par_HRUtype::KS);
    prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::LINL_RES:
    BaseOut = (prev_Grou * get_par(par_HRUtype::L)) * get_par(par_HRUtype::KS) ;
    prev_Grou = (prev_Grou * get_par(par_HRUtype::L)) + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::LINBY_RES:
    BaseOut = prev_Grou * get_par(par_HRUtype::KS) + get_par(par_HRUtype::D_BYPASS) * ((1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC));
    prev_Grou = prev_Grou  + (1 - get_par(par_HRUtype::D_BYPASS)) * (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - prev_Grou * get_par(par_HRUtype::KS);

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::POW_RES:
      BaseOut = std::pow(prev_Grou, get_par(par_HRUtype::B_EXP)) * get_par(par_HRUtype::KS);
      BaseOut = std::min(BaseOut, prev_Grou);
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

      set_varValue(BaseOut, tstRM, ts_type::BASF);
      set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::EXP_RES:
    if(get_par(par_HRUtype::B_EXP) != 0) {
      BaseOut = get_par(par_HRUtype::KS)*std::exp(prev_Grou / get_par(par_HRUtype::B_EXP));
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

      set_varValue(BaseOut, tstRM, ts_type::BASF);
      set_varValue(prev_Grou, tstRM,ts_type::GROS);
    } else{
      BaseOut = 0.0;
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;
      set_varValue(BaseOut, tstRM, ts_type::BASF);
      set_varValue(prev_Grou, tstRM,ts_type::GROS);
    }
    break;

  case gs_STORtype::LIN_2SE:
    BaseOut = prev_GroS2 * get_par(par_HRUtype::KS2);
    BaseOut_1 = prev_GroS1 * get_par(par_HRUtype::KS);

    prev_GroS1 = prev_GroS1 + ((1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC)) - BaseOut_1;
    prev_GroS2 = prev_GroS2 + BaseOut_1 - BaseOut;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    // set_varValue(prev_Grou1, tstRM,ts_type::GROS);
    // set_varValue(prev_Grou2, tstRM,ts_type::GROS);
    set_varValue(prev_GroS1 + prev_GroS2, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::LIN_2PA:
    BaseOut_1 = prev_GroS1 * get_par(par_HRUtype::KS);
    prev_GroS1 = prev_GroS1 + (1 - get_par(par_HRUtype::ADIV)) * get_par(par_HRUtype::ALPHA) * get_dta(tstRM, ts_type::PERC) - BaseOut_1;

    BaseOut_2 = prev_GroS2 * get_par(par_HRUtype::KS2);
    prev_GroS2 = prev_GroS2 + (1 - get_par(par_HRUtype::ADIV)) * (1 - get_par(par_HRUtype::ALPHA)) * get_dta(tstRM, ts_type::PERC) - BaseOut_2;

    BaseOut = BaseOut_1 + BaseOut_2;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    // set_varValue(prev_Grou1, tstRM,ts_type::GROS);
    // set_varValue(prev_Grou2, tstRM,ts_type::GROS);
    set_varValue(prev_GroS1 + prev_GroS2, tstRM,ts_type::GROS);

    break;

  case gs_STORtype::FLEX_RES:
   if(get_par(par_HRUtype::THR) >= prev_Grou) {
      //lower outlet working
      BaseOut_1 = prev_Grou * get_par(par_HRUtype::KS);
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut_1;
      set_varValue(BaseOut_1, tstRM, ts_type::BASF);

    } else {
    //   //lower and upper outlets working
    //
    //   //upper outlet
      // BaseOut_1=0.0;
      BaseOut_2 = get_par(par_HRUtype::KS2) * (prev_Grou - get_par(par_HRUtype::THR));
      prev_Grou = prev_Grou - BaseOut_2;
    //
    //   //lower outlet
    // BaseOut_2=0.0;
      BaseOut_1 = prev_Grou * get_par(par_HRUtype::KS);
      prev_Grou = prev_Grou - BaseOut_1;
    //
    //
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC);
    //
      BaseOut = BaseOut_1 + BaseOut_2;
    //
      set_varValue(BaseOut, tstRM, ts_type::BASF);
    }

    set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::EXP_LOG:
    if(get_par(par_HRUtype::B_EXP) != 0) {
      BaseOut = get_par(par_HRUtype::B_EXP) * std::log(1+std::exp(prev_Grou / get_par(par_HRUtype::B_EXP)));
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

    } else{
      BaseOut = 0.0;
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;
    }

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou, tstRM,ts_type::GROS);

    break;
  }



  return;
}

/** \brief Updates the series of fast response described by linear reservoirs
 *
 */

void single_HMunit::fast_response() {

  numberSel helpFastOut = 0.0, help_State =0.0;

  //  help_State = get_stateFastres(0);
// std::cout << " help_nmbr before" << help_nmbFR << std::endl;
  for(ifrb=0; ifrb<help_nmbFR; ifrb++) {
    help_State = get_stateFastres(ifrb);
    helpFastOut = get_par(par_HRUtype::KF) * help_State;
    // help_State = help_State - helpFastOut;
    if(ifrb == 0) {
      help_State = help_State + get_par(par_HRUtype::ADIV) * get_dta(tstRM,ts_type::PERC)- helpFastOut;
    } else help_State = help_State + get_outFastRes((ifrb-1))- helpFastOut;
    set_stateFastRes(help_State,ifrb);
    set_outFastRes(helpFastOut,ifrb);
    // std::cout << ifrb << " ifrb " << get_stateFastres(ifrb) << " state " << ifrb << " ifrb " << helpFastOut << " adivresc " << get_par(par_HRUtype::ADIV) * get_dta(tstRM,ts_type::PERC) << std::endl;
  }

  set_varValue(helpFastOut,tstRM,ts_type::DIRR);
  //  set_varValue(help_State,tstRM,ts_type::SURS);
  return ;
}

/** \brief updates states in interception storage in single pdm unit
 *
 *  updating of states in canopy and stem-trunk storage in the single pdm unit
 *  calculating related fluxes
 *
 */
void single_HMunit::interception_NoSnow(interception_STORtype _intrc_STORAGE) {


  //
  switch(_intrc_STORAGE) {

  case interception_STORtype::Rutter_Gash:

  numberSel CanOut = 0.0, StemOut = 0.0, OverflowCan = 0.0, OverflowStem, EvapCanop = 0.0, EvapStem = 0.0, Througf = 0.0;

  OverflowCan = std::max((prevCanS - get_par(par_HRUtype::CAN_ST)),0.0);
  //!< VIC model for canopy evaporation (prevCanS/ get_par(par_HRUtype::CAN_ST))^(2/3)
  prevCanS = prevCanS - OverflowCan;
  EvapCanop = std::min(std::pow(((prevCanS) / get_par(par_HRUtype::CAN_ST)),2/3),prevCanS);

  numberSel help_EvapCanop = update_ETDEMAND(EvapCanop, false);
  et_demand = update_ETDEMAND(EvapCanop, true);
  EvapCanop = help_EvapCanop;


  //
  //
  // if(EvapCanop >= et_demand) {
  //   EvapCanop = et_demand;
  //   et_demand = 0.0;
  // }
  // else {
  //   et_demand = et_demand - EvapCanop;
  // }

  prevCanS = prevCanS - EvapCanop;

  CanOut = std::min((prevCanS / get_par(par_HRUtype::CAN_ST) * EvapCanop),prevCanS);
  prevCanS = prevCanS - CanOut;
  set_varValue(prevCanS, tstRM, ts_type::CANS);

  prevCanS =  prevCanS + get_par(par_HRUtype::CDIV) * (get_dta(tstRM, ts_type::PREC) + get_par(par_HRUtype::CDIV) *get_dta(tstRM, ts_type::MELT));

  //!< VIC model for canopy evaporation (prevCanS/ get_par(par_HRUtype::CAN_ST))^(2/3)
  OverflowStem = std::max((prevSteS - get_par(par_HRUtype::STEM_ST)),0.0);
  prevSteS = prevSteS - OverflowStem;

  EvapStem = std::min(std::pow(((prevSteS) / get_par(par_HRUtype::STEM_ST)),(2/3)), prevSteS);

  numberSel help_EvapStem = update_ETDEMAND(EvapStem, false);
  et_demand = update_ETDEMAND(EvapStem, true);
  EvapStem = help_EvapStem;

  // if(EvapStem >= et_demand) {
  //   EvapStem = et_demand;
  //   et_demand = 0.0;
  // }
  // else {
  //   et_demand = et_demand - EvapStem;
  // }

  prevSteS = prevSteS - EvapStem;

  StemOut = std::min((prevCanS) / get_par(par_HRUtype::CAN_ST) * EvapStem, prevSteS);
  prevSteS = prevSteS - StemOut;
  set_varValue(prevSteS, tstRM, ts_type::STES);

  prevSteS = prevSteS + get_par(par_HRUtype::SDIV) * (get_dta(tstRM, ts_type::PREC) + get_par(par_HRUtype::SDIV) *get_dta(tstRM, ts_type::MELT)) + (1 - get_par(par_HRUtype::CSDIV)) * (CanOut + OverflowCan);

  set_varValue((CanOut + OverflowCan), tstRM, ts_type::CANF);
  set_varValue(EvapCanop, tstRM, ts_type::EVAC);

  set_varValue((StemOut + OverflowStem), tstRM, ts_type::STEF);
  set_varValue(EvapStem, tstRM, ts_type::EVAS);

  Througf = (OverflowCan + CanOut) * get_par(par_HRUtype::CSDIV) + StemOut + OverflowStem;

  set_varValue(Througf, tstRM, ts_type::TROF);

  set_varValue((get_dta(tstRM, ts_type::CANS) + get_dta(tstRM, ts_type::STES)),tstRM,ts_type::INTS);


  break ;
  }

  return ;
}



/** \brief updates states in interception storage in single HRU unit
 *
 *  updating of states in canopy and stem-trunk storage in the single hru unit
 *  calculating related fluxes
 *
 */
void single_HMunit::interception_WithSnow(interception_STORtype _intrc_STORAGE) {

  switch(_intrc_STORAGE) {

  case interception_STORtype::Rutter_Gash:
  //  numberSel CanOut = 0.0, StemOut = 0.0, OverflowCan = 0.0, OverflowStem, EvapCanop = 0.0, EvapStem = 0.0, Througf = 0.0;
  numberSel OverflowCan = 0.0, OverflowStem= 0.0, EvapCanop = 0.0, EvapStem = 0.0, Througf = 0.0;

  // OverflowCan = std::max((prevCanS - get_par(par_HRUtype::CAN_ST)),0.0);
  // prevCanS = prevCanS - OverflowCan;
  // prevCanS = prevCanS - OverflowCan;
  //!< VIC model for canopy evaporation (prevCanS/ get_par(par_HRUtype::CAN_ST))^(2/3)
  //the sublimation of snow
  EvapCanop = std::min(pow(((prevCanS) / get_par(par_HRUtype::CAN_ST)),2/3),prevCanS);

  numberSel help_EvapCanop = update_ETDEMAND(EvapCanop, false);
  et_demand = update_ETDEMAND(EvapCanop, true);
  EvapCanop = help_EvapCanop;

  // if(EvapCanop >= et_demand) {
  //   EvapCanop = et_demand;
  //   et_demand = 0.0;
  //   }
  //  else {
  //    et_demand = et_demand - EvapCanop;
  //    }

  //EvapCanop = 0.0;
  prevCanS = prevCanS - EvapCanop;
  //  CanOut = (prevCanS - OverflowCan) / get_par(par_HRUtype::CAN_ST) * EvapCanop;
  //  prevCanS =  prevCanS + (get_dta(tstRM, ts_type::SNOW) + get_dta(tstRM, ts_type::MELT))- OverflowCan - CanOut - EvapCanop;
  set_varValue(prevCanS, tstRM, ts_type::CANS);

  prevCanS =  prevCanS + get_par(par_HRUtype::CDIV) * get_dta(tstRM, ts_type::MELT);

  OverflowStem = std::max((prevSteS - get_par(par_HRUtype::STEM_ST)),0.0);
  prevSteS = prevSteS - OverflowStem;
  //!< VIC model for canopy evaporation (prevCanS/ get_par(par_HRUtype::CAN_ST))^(2/3)
  EvapStem = std::min(pow(((prevSteS) / get_par(par_HRUtype::STEM_ST)),(2/3)), prevSteS);

  numberSel help_EvapStem = update_ETDEMAND(EvapStem, false);
  et_demand = update_ETDEMAND(EvapStem, true);
  EvapStem = help_EvapStem;
  //
  // if(EvapStem >= et_demand) {
  //   EvapStem = et_demand;
  //   et_demand = 0.0;
  // }
  // else {
  //   et_demand = et_demand - EvapStem;
  // }
  //

  //EvapStem = 0.0;
  prevSteS = prevSteS - EvapStem;
  //  StemOut = (prevCanS - OverflowCan) / get_par(par_HRUtype::CAN_ST) * EvapStem;
  //  prevSteS = prevSteS + get_par(par_HRUtype::SDIV) * (get_dta(tstRM, ts_type::PREC) + get_dta(tstRM, ts_type::MELT)) + (1 - get_par(par_HRUtype::CSDIV)) * (CanOut + OverflowCan) - EvapStem - StemOut;
  set_varValue(prevSteS, tstRM, ts_type::STES);

  prevSteS = prevSteS + get_par(par_HRUtype::SDIV) * get_dta(tstRM, ts_type::MELT);

  //  set_varValue((CanOut + OverflowCan), tstRM, ts_type::CANF);
  set_varValue(0.0, tstRM, ts_type::CANF);

  set_varValue(EvapCanop, tstRM, ts_type::EVAC);

  set_varValue((OverflowCan + OverflowStem), tstRM, ts_type::STEF);
  // set_varValue(0.0, tstRM, ts_type::STEF);
  set_varValue(EvapStem, tstRM, ts_type::EVAS);

  //  Througf = (OverflowCan + CanOut) * get_par(par_HRUtype::CSDIV) + StemOut + OverflowStem;

  set_varValue(Througf, tstRM, ts_type::TROF);
  // set_varValue((prevCanS + prevSteS),tstRM,ts_type::INTS);
  set_varValue((get_dta(tstRM, ts_type::CANS) + get_dta(tstRM, ts_type::STES)),tstRM,ts_type::INTS);

  //   std::cout <<  " prevCanS " << prevCanS << " EvapCanop " << EvapCanop << "\n";
  break;
  }

  return ;
}


/** \brief Updates interception and snow buffer
 *
 */
void single_HMunit::snow_melt(){

  numberSel Snow_melt = 0.0, Snoww = 0.0;

  Snoww = get_dta(tstRM, ts_type::SNOW);

  if(get_dta(tstRM, ts_type::TEMP) > get_par(par_HRUtype::TMEL)) {
    Snow_melt = std::min(get_par(par_HRUtype::DDFA) * (get_dta(tstRM, ts_type::TEMP) - get_par(par_HRUtype::TMEL)), prevSnoS);
    // if(Snow_melt <0) std::cout  <<"negative melt  " << Snow_melt << "  "<< prevSnoS << " " << get_par(par_HRUtype::TMEL) << " " <<get_dta(tstRM, ts_type::TEMP) <<" \n";
  } else Snow_melt = 0.0;

  if ((prevSnoS + Snoww - Snow_melt)<0){
    prevSnoS = 0.0;
    Snow_melt = prevSnoS + Snoww;
  } else prevSnoS = prevSnoS + Snoww - Snow_melt;

  //   std::cout << " prevSnoS " << prevSnoS << " MELT " << Snow_melt<< " ";


  // set_varValue(prevSnoS,tstRM,ts_type::SNOW);
  set_varValue(Snow_melt,tstRM,ts_type::MELT);


  return ;
}

/** \brief Updates interception and snow buffer
 *
 */
void single_HMunit::interception_snow() {

  numberSel Snow_melt = 0.0, Snoww = 0.0;

  if(get_dta(tstRM, ts_type::TEMP) < get_par(par_HRUtype::TETR)) {
    Snoww = get_dta(tstRM, ts_type::PREC);
    set_varValue(Snoww,tstRM,ts_type::SNOW);
    snow_melt();
    interception_WithSnow(intrc_STORAGE);
    //    std::cout << " snow " << Snoww << " \n";
  } else {
    Snoww = 0.0;
    set_varValue(Snoww,tstRM,ts_type::SNOW);
    snow_melt();
    interception_NoSnow(intrc_STORAGE);
  }

  // set_varValue(Snoww,tstRM,ts_type::SNOW);
//
//   if(get_dta(tstRM, ts_type::TEMP) > get_par(par_HRUtype::TMEL)) {
//     Snow_melt = std::min(get_par(par_HRUtype::DDFA) * (get_dta(tstRM, ts_type::TEMP) - get_par(par_HRUtype::TMEL)), prevSnoS);
//     // if(Snow_melt <0) std::cout  <<"negative melt  " << Snow_melt << "  "<< prevSnoS << " " << get_par(par_HRUtype::TMEL) << " " <<get_dta(tstRM, ts_type::TEMP) <<" \n";
//   } else Snow_melt = 0.0;
//
//   if ((prevSnoS + Snoww - Snow_melt)<0){
//     prevSnoS = 0.0;
//     Snow_melt = prevSnoS + Snoww;
//   } else prevSnoS = prevSnoS + Snoww - Snow_melt;
//
//   //   std::cout << " prevSnoS " << prevSnoS << " MELT " << Snow_melt<< " ";
//
//
//   // set_varValue(prevSnoS,tstRM,ts_type::SNOW);
//   set_varValue(Snow_melt,tstRM,ts_type::MELT);


  return ;
}

/** \brief Updates fluxes and states in for loop for all time intervals, set zeros to internal help state and initial states
 *
 */
void single_HMunit::run_HB() {

  unsigned Numdta;

  Numdta =  get_numdta();
  // if( !(Numdta >0)) {
  //   std::cout << std::endl << "There is an error in data loadings." << std::endl;
  //   std::cout << "It is impossible to calculate the basic for loop in runHB function" << std:: endl << "It is controlled by the length " << Numdta << "." << std::endl;
  //   std::exit(EXIT_FAILURE);
  // }

  // std::cout << "blalal " << help_nmbFR <<"\n";
  numberSel helprm=0.0;
  for(tstRM=0; tstRM < Numdta ; tstRM++) {
    et_demand = get_dta(tstRM,ts_type::PET);
    // std::cout << " beg "<< et_demand << "\n";
    interception_snow();//
    // std::cout << " interception "<< et_demand << " evac " << get_dta(tstRM, ts_type::EVAC) << " evas " << get_dta(tstRM, ts_type::EVAS) <<"\n";
    surface_retention();//
    // std::cout << " surf ret "<< et_demand << " ewsr " << get_dta(tstRM,ts_type::ETSW) <<"\n";
    // std::cout << tstRM << "\n\n";
    soil_buffer(soil_STORAGE);//
    slow_response(gs_STORAGE);
    fast_response();
    helprm = (get_dta(tstRM,ts_type::BASF) + get_dta(tstRM,ts_type::DIRR));
    set_varValue(helprm ,tstRM,ts_type::TOTR);
    upadate_actualET();

    //    std::cout <<(get_dta(tstRM,ts_type::BASF) + get_dta(tstRM,ts_type::DIRR)) << " "<< get_dta(tstRM,ts_type::BASF) << " "<< get_dta(tstRM,ts_type::DIRR)<< "\n";
  }
  //  std::cout << "prev_Ground storage before zeros " << prev_Grou << std::endl;
  //  std::cout << "prevCanS  " << prevCanS << std::endl;
  //  std::cout << "prevSteS " << prevSteS << std::endl;
  //  std::cout << "prevSnoS " << prevSnoS << std::endl;
  //  std::cout << "prev_SurS " << prev_SurS << std::endl;
  //  std::cout << "prev_Soil " << prev_Soil << std::endl;
  //  std::cout << "prev_Grou " << prev_Grou << std::endl;
  //  std::cout << "precolation " << get_dta((tstRM-3),ts_type::PREC) << std::endl;
  //  std::cout << "Tot RM " << get_dta((tstRM-1),ts_type::TOTR) << std::endl;
  //  std::cout << "Tot Soil " << get_dta((tstRM-3),ts_type::SOIS) << std::endl;

  //  set_ZeroinitStates(help_nmbFR);
  //  std::cout << "prev_Grouendfrb  " << prev_Grou << std::endl;
  set_ZeroStates();
  //  std::cout << "prev_Ground storage after zeros " << prev_Grou << std::endl;

  return ;

}

/** \brief initialization of all hdata of all state variables and fluxes using the value val
 *
 * \param value of initialization
 *
 */
void single_HMunit::init_inputs(numberSel val, unsigned numDTA) {

  hdata dta(val, numDTA);
    // std::cout << "Initializing the hdata setup " << get_numdta() << " ups numDTA "<< numDTA << std::endl;
  set_data(dta,ts_type::PREC);
  // std::cout << " prec ok\n";
  set_data(dta,ts_type::AET);
  set_data(dta,ts_type::PET);
  //!< Snow melting process DDM
  set_data(dta,ts_type::TEMP);
  set_data(dta,ts_type::SNOW);
  set_data(dta,ts_type::MELT);
  //!< Interception variables, evapotranspiration
  set_data(dta,ts_type::CANF);
  set_data(dta,ts_type::CANS);
  set_data(dta,ts_type::EVAC);
  set_data(dta,ts_type::EVAS);
  set_data(dta,ts_type::STEF);
  set_data(dta,ts_type::STES);
  set_data(dta,ts_type::TROF);
  set_data(dta,ts_type::INTS);
  //!< Surface retention
  set_data(dta,ts_type::SURS);
  set_data(dta,ts_type::ETSW);
  set_data(dta,ts_type::PREF);
  // std::cout << " precf ok\n";
  //!< soil percolation
  set_data(dta,ts_type::EVBS);
  set_data(dta,ts_type::SOIS);
  set_data(dta,ts_type::PERC);
  //!< Groundwater variables
  set_data(dta,ts_type::BASF);
  set_data(dta,ts_type::GROS);
  //!< Fast runoff response
  set_data(dta,ts_type::DIRR);
  set_data(dta,ts_type::TOTR);
  // std::cout << " tor ok\n";

  // get_numdta();

     // std::cout << "The total number of initialized time intervals is " << get_numdta() << " ." << std::endl;

  return ;

}

/** \brief Loading of precipitatin and temperature data
 *
 * \param precipitation data same size as temperature data
 * \param temperature data same size as precipitation data
 *
 */
void single_HMunit::load_data_PT(const hdata& prec_input, const hdata& temp_input, const numberSel& val,const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay) {

  unsigned helpnumDTA;

  helpnumDTA = prec_input.size();

  // std::cout << "ups size " << prec_input.size() << "\n";

  // if(helpnumDTA != temp_input.size()) {
  //   std::cout << "Different number of time intervals in precipitation input " << prec_input.size() \
  //             << " the temperature input has " << temp_input.size() << " inputs." << std::endl;
  //   std::exit(EXIT_FAILURE);
  //
  // }
  //  numberSel val = -99999.9;
  //  numberSel val = 0;
  // std::cout << "\n init inputs 1\n";
  init_inputs(val, helpnumDTA);
  // std::cout << "\n init inputs 2\n";
  set_data_prec_temp(prec_input,temp_input);
  // std::cout << "\n PT data inputs 1\n";
  set_calender(inYear,inMonth,inDay,helpnumDTA);
  // std::cout << "\n class 1\n";

  return ;

}

/** \brief Setting of initial date and than set year month day a Jday for each time interval
 *
 * \param year of the first date
 * \param month of the first date
 * \param day of the first date
 * \param number of time intervals
 *
 */
void single_HMunit::set_calender(const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay,const unsigned& initNumTS) {

  hyd_dta.s_initDate(inYear,inMonth,inDay, initNumTS);
  hyd_dta.s_calender();
  //  hyd_dta.p_calender();

  return ;

}

/** \brief Setting the PET Varibales
 *
 * \param latitude
 * \param ID of PEt type
 *
 */
void single_HMunit::set_PetVars(const numberSel& newLatitude, const pet_Type& newPeType) {

  hyd_dta.s_Pet_Pars( newLatitude, newPeType);

  return ;
}

/** \brief Calculates the PET default OUDIN methods is used
 *
 */

void single_HMunit::calc_Pet() {

  hyd_dta.calc_Pet();

  return ;

}

/** \brief Loading the new param values
 *
 * \param vector of pair values and ID of parameter
 *
 */

void single_HMunit::set_paramsToSim(std::vector<std::pair<numberSel,par_HRUtype>>& parsToLoad) {

  //    std::cout << std::endl << "Params before loadings:" << std::endl;
  //    par_HRU.p_param();
  par_HRU.s_parLoadToCalib(parsToLoad);
  //    std::cout << std::endl << "Params after loadings:" << std::endl;
  //    par_HRU.p_param();

  return;

}

/** \brief Set default values of parameters and
 *
 *  \param if true print values to cout
 *
 */
void single_HMunit::p_defaultParams(bool Bprint) {

  par_HRU.s_default();
  if(Bprint)  par_HRU.p_param();

  return ;

}

/** \brief Set zero on all state variables
 *
 *
 */
void single_HMunit::set_ZeroStates() {

  prev_Soil = 0.0;
  prev_Grou = 0.0;
  prevCanS = 0.0;
  prevSteS = 0.0;
  prevSnoS = 0.0;
  prev_SurS = 0.0;
  prev_GroS1 = 0.0;
  prev_GroS2 = 0.0;

  help_nmbFR = get_nmbFastRes();
  //  std::cout << "\nFast runoff response has " << help_nmbFR << " reservoirs." << std::endl;
  set_ZeroinitStates(help_nmbFR);

  return ;

}

/** \brief Prints hyd_dta into the file
 *
 * \param name of file with path in string to save data in
 *
 */
void single_HMunit::print_OutputToFile(const std::string& Filet) {

  hyd_dta.printDataToFile(Filet);

  //  std::ofstream outfilet;
  //  outfilet.open (Filet.c_str());
  //
  //  if(outfilet.is_open()) {
  //    outfilet.precision(6);
  //
  //    unsigned Numdta;
  //    Numdta =  get_numdta();
  //    outfilet << std::right;
  //    char setfiller = ' ';
  //    for(unsigned Ts=0; Ts<Numdta; Ts++) {
  //      for(const auto& it : all_caDT) {
  //        outfilet << std::setprecision(0);
  //        outfilet << std::setw(2);
  //        outfilet << hyd_dta.g_calDta(it,Ts) << " ";
  //      }
  //
  //      for(const auto& it : all_ts) {
  //        outfilet << std::setprecision(2);
  //        outfilet << std::fixed << std::setw(7);
  //        outfilet << std::setfill(setfiller) << std::right;
  //        outfilet << get_dta(Ts, it) << " ";
  //      }
  //      outfilet << std::endl;
  //    }
  //    outfilet.close();
  //
  //  } else {
  //
  //    std::cout << "There is a error in opening and writing to the file with path and name: " << Filet.c_str() << std::endl;
  //
  //  }

  return ;

}


/** \brief Read data for vallarray of input vars from the file
 *
 * \param name of file with path in string to read data from
 *
 */
void single_HMunit::read_InputFromFile(const std::string& Filet) {

  std::ifstream infilet;
  numberSel helpInitYear = 0.0, helpInitMonth = 0.0,helpInitDay = 0.0;
  unsigned numLinesInFile = 0;
  std::string line;

  hdata helpTemp, helpPrec;
// std::cout << "strating opening\n";
  infilet.open(Filet.c_str());
  if(infilet.is_open()) {
    while (std::getline(infilet, line)) ++numLinesInFile;
     // std::cout << "Number of lines in file: " << Filet.c_str() << " equals to " << numLinesInFile << "." << std::endl;
    if(line.empty()) numLinesInFile--;// if last line is empty.
        // std::cout << "Number of lines in file: " << Filet.c_str() << " equals to " << numLinesInFile << "." << std::endl;
    infilet.clear();
    infilet.seekg (0, std::ios::beg);
    infilet >> helpInitYear >> helpInitMonth >>helpInitDay;
        // std::cout << helpInitYear << " " << helpInitMonth << " " << helpInitDay << std::endl;
    numLinesInFile--;
    helpTemp.resize(numLinesInFile);
    helpPrec.resize(numLinesInFile);
        // std::cout << "foor loop Number of lines in file: " << Filet.c_str() << " equals to " << numLinesInFile << "." << std::endl;
    for(unsigned it=0; it<numLinesInFile; it++) {
      std::string lineVal;
      if(it==0) getline(infilet,lineVal);// reading the first line with year month  day
      getline(infilet, lineVal);
      std::stringstream valueBuffer(lineVal);
      valueBuffer >> helpTemp[it] >> helpPrec[it];
    }
        // std::cout <<std::endl << helpPrec[0] << " " << helpTemp[0] << " size " << helpPrec.size() << std::endl;
    infilet.close();
// std::cout << "load_data";
    load_data_PT(helpPrec,helpTemp,0.0,helpInitYear,helpInitMonth,helpInitDay);

    // std::cout << "\nload_data loade";

  }
  // else {
  //
  //   std::cout << "There is a error in opening and reading from the file with path and name: " << Filet.c_str() << std::endl;
  //
  // }

  return ;

}

/** \brief Setter for area of HRU
 *
 * \param area size of HRU in m2
 *
 */
void single_HMunit::set_Area(numberSel area) {

  Area = area;

  return ;

}

/** \brief Getter for area of HRU
 *
 * \param area size of HRU in m2
 *
 */
numberSel single_HMunit::get_Area(){

  return Area;

}

/** \brief Getter for single time series of HRU
 *
 * \param _tsType type of time serie of HRU
 *
 */
hdata single_HMunit::getSingleHruTsDta(const ts_type& _tsType){

 return hyd_dta.get_HbTsData(_tsType);

}

/** \brief Getter for all time series of HRU
 *
 *
 */
data_HB_1d single_HMunit::getAllData(){

 return hyd_dta;

}

/** \brief Setter for ID of HRU
 *
 * \param IdToSet IS string of HRU
 *
 */
void single_HMunit::setIdHru(const std::string& IdToSet) {

  IdHru = IdToSet;

 return ;
}

/** \brief Getter for ID of HRU
 *
 *
 */
std::string single_HMunit::getIdHru(){

 return (IdHru);

}

/** \brief Setter for calendar dat of HRU
 *
 * \param yyear the data of years for caldata
 * \param mmonth the data of mmonths for caldata
 * \param dday the data of days for caldata
 *
 */
void single_HMunit::load_calData(const caldata& yyear, const caldata& mmonth, const caldata& dday) {

  hyd_dta.loadCalData(yyear, mmonth, dday);

  return ;

}

/** \brief Getter for number of parameters of HRU
 *
 *
 */
unsigned single_HMunit::get_numPars(){

 return (par_HRU.g_numPars());

}

/** \brief Prints all names and values of parameters of given  HRU
 *
 *
 */
void single_HMunit::print_Pars() {

  // std::cout << "\nThe HRU ID " << getIdHru() << std::endl;

  par_HRU.p_param();

  return ;

}

void single_HMunit::set_GStype(gs_STORtype _gs_STORtype) {

  gs_STORAGE = _gs_STORtype;

  return ;

}

gs_STORtype single_HMunit::get_GStype() {
  return gs_STORAGE;

}

numberSel single_HMunit::update_ETDEMAND(const numberSel& ET, bool ET_demand){

  numberSel et_corrected = 0.0, et_demand_updated = 0.0, value = 0.0;
  // if(std::isnan(et_demand)) et_demand = 0.0;

  if(ET >= et_demand) {
    et_corrected = et_demand;
    et_demand_updated = 0.0;
  }
  else {
    et_demand_updated = et_demand - ET;
    et_corrected = ET;
  }

  if(ET_demand) {
    value = et_demand_updated;
  } else {
      value  = et_corrected;
         }

  return value;
}

// void single_HMunit::print_GStype() {
//
//   switch(gs_STORAGE) {
//
//   case gs_STORtype::LIN_RES:
//
//     std::cout << "The gs_STORE is a LIN reservoir." << std::endl;
//
//     break;
//
//   case gs_STORtype::LINL_RES:
//
//     std::cout << "The gs_STORE is a LINL reservoir." << std::endl;
//
//     break;
//
//   case gs_STORtype::LINBY_RES:
//
//     std::cout << "The gs_STORE is a LINBY reservoir." << std::endl;
//
//     break;
//
//   case gs_STORtype::LIN_2SE:
//
//     std::cout << "The gs_STORE is a LIN_2SE reservoir." << std::endl;
//
//     break;
//
//   case gs_STORtype::LIN_2PA:
//
//     std::cout << "The gs_STORE is a LIN_2PA reservoir." << std::endl;
//
//     break;
//
//   case gs_STORtype::POW_RES:
//
//     std::cout << "The gs_STORE is a POW reservoir." << std::endl;
//
//     break;
//
//   case gs_STORtype::EXP_RES:
//
//     std::cout << "The gs_STORE is a EXP reservoir." << std::endl;
//
//     break;
//
//   case gs_STORtype::FLEX_RES:
//
//     std::cout << "The gs_STORE is a FLEX reservoir." << std::endl;
//
//     break;
//
//   }
//
// }

void single_HMunit::set_soilStorType(soil_STORtype _soil_STORtype) {

  soil_STORAGE = _soil_STORtype;
  //print_soilStorType();

  return ;

}

soil_STORtype single_HMunit::get_soilStorType() {
  return soil_STORAGE;

}

void single_HMunit::print_soilStorType() {

  switch(soil_STORAGE) {

  case soil_STORtype::PDM:
    std::cout << "The soil_STOR is a PDM reservoir." << std::endl;
    break;

  case soil_STORtype::COLLIE_V2:
    std::cout << "The soil_STOR is a COLLIE_V2 reservoir." << std::endl;
    break;

  case soil_STORtype::NEW_ZEALAND:
    std::cout << "The soil_STOR is a NEW_ZEALAND reservoir." << std::endl;
    break;

  case soil_STORtype::GR4J:
    std::cout << "The soil_STOR is a GR4J reservoir" << std::endl;
    break;

  case soil_STORtype::SBROOK_V1:
    std::cout << "The soil_STOR is a SBROOK_V1 reservoir" << std::endl;
    break;

  case soil_STORtype::HILLSLOPE:
    std::cout << "The soil_STOR is a HILLSLOPE reservoir" << std::endl;
    break;

  case soil_STORtype::PLATEAU:
    std::cout << "The soil_STOR is a PLATEAU reservoir" << std::endl;
    break;

  case soil_STORtype::PDM2:
    std::cout << "The soil_STOR is a PDM2 reservoir." << std::endl;
    break;

  }

}

std::vector<numberSel> single_HMunit::water_balance(numberSel next_soil, numberSel val, std::vector<numberSel> vals) {


  for(unsigned i=0; i<vals.size(); i++) {
    if(next_soil >= 0) {
      if(next_soil < vals[i]) {
        vals[i] = next_soil;
      }
      next_soil = next_soil - vals[i];
      val = val + vals[i];
    }

  }

  std::vector<numberSel> updated_vals{next_soil, val};

  return updated_vals;

}

void single_HMunit::upadate_actualET() {

  numberSel aet = 0.0;

  aet = get_dta(tstRM, ts_type::EVAC) + get_dta(tstRM, ts_type::EVAS) + get_dta(tstRM, ts_type::ETSW) + get_dta(tstRM, ts_type::EVBS);
  set_varValue(aet, tstRM,ts_type::AET);

}


void single_HMunit::set_inteceptionType(interception_STORtype _intrc_STORAGE){

  intrc_STORAGE = _intrc_STORAGE;

}

interception_STORtype single_HMunit::get_intercetionStorType(){

  return intrc_STORAGE;

}

