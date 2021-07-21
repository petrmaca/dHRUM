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
  help_nmbFR(0),
  ifrb(0),
  Area(0),
  IdHru() {

  set_nmbFastres(10);
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
  p_defaultParams(false);

//  std::cout << "INITprevDR " << prev_Grou << std::endl;

  tstRM = 0;
  gs_STORAGE = gs_STORtype::LIN_2SE;

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
help_nmbFR(0),
ifrb(0),
Area(0),
IdHru() {

  tstRM = other.tstRM;//!< The counter for main loop in run model
  par_HRU = other.par_HRU;//!< The parameters in PDM instances
  hyd_dta = other.hyd_dta;//!< The data of all time series of hydrological variables
  prev_Soil = other.prev_Soil;//!< The helper variable for updating soil storage
  prev_Grou = other.prev_Grou;//!< The helper variable for updating groundwater storage
  prevCanS = other.prevCanS;//!<  The helper variable for Canopy interception storage
  prevSteS = other.prevSteS;//!<  The helper variable for Stem interception storage
  prevSnoS = other.prevSnoS;//!<  The helper variable for Snow storage
  prev_SurS = other.prev_SurS;//!< The helper variable for updating surface storage
  help_nmbFR = other.help_nmbFR;//!< The helper for number of fast reservoirs
  ifrb = other.ifrb;//!< For loop counter
  Area = other.Area;//!< The area of HM unit in m2
  IdHru = other.IdHru;

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
    help_nmbFR = rhs.help_nmbFR;//!< The helper for number of fast reservoirs
    ifrb = rhs.ifrb;//!< For loop counter
    Area = rhs.Area;//!< The area of HM unit in m2
    IdHru = rhs.IdHru;

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
  if (get_dta(tstRM, ts_type::TEMP) > 0.0) {
    EvapSR = 0.0824 * std::pow(get_dta(tstRM, ts_type::TEMP),1.289);
  } else EvapSR = 0.0;

  if (EvapSR > prev_SurS) {
    EvapSR = prev_SurS;
  }

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

  set_varValue(EvapSR, tstRM, ts_type::AET);
  set_varValue(RetOut,tstRM,ts_type::PREF);

  return ;
}


/** \brief updates states in soil buffer in single pdm unit
 *
 *  updating of states in soil buffer in single pdm unit
 *  calculating related fluxes
 *
 */
void single_HMunit::soil_buffer() {
  //     if(tstRM == 0){
  //     if(prev_Soil > get_par(par_HRUtype::SMAX)) {
  //     // overFl0 = prev_Soil - get_par(par_HRUtype::SMAX);
  //     prev_Soil = get_par(par_HRUtype::SMAX);
  //   }
  // }

  numberSel c_init = 0.0, overFl1 = 0.0, ppInf = 0.0, overFl2 = 0.0, c_prop = 0.0, next_soil = 0.0, overFL = 0.0, evap = 0.0, aet = 0.0, c_contr=0.0;
  ////  overflow from previous day
  // if(prev_Soil > get_par(par_HRUtype::SMAX)) {
  //   overFl0 = prev_Soil - get_par(par_HRUtype::SMAX);
  //   prev_Soil = get_par(par_HRUtype::SMAX);
  // }
  // Estimation of Soil Water Depth C using total soil basin Storage S from previous day see Moore description of PDM HESS 2007,
  //  Wood and bascics from 1992
  // Eric F. Wood, D. P. Lettenmaier, V. G. Zartarian A land-surface hydrology parameterization with subgrid variability for general circulation models
//   // eq 3a or 18a Vic paper
// // //// to test
//     if((get_dta(tstRM, ts_type::PREF) + prev_Soil)>=get_par(par_HRUtype::SMAX)) {
//      overFl1 = get_dta(tstRM, ts_type::PREF) - get_par(par_HRUtype::SMAX) + prev_Soil;
//      overFl2 =0.0;
//      next_soil = 0.0;
//      // std::cout <<"\nupod";
//      evap =   std::min(static_cast<numberSel>(get_dta(tstRM, ts_type::PREF) + prev_Soil - overFl1), static_cast<numberSel> (get_dta(tstRM, ts_type::PET) * (1 - pow(((get_par(par_HRUtype::SMAX)) / get_par(par_HRUtype::SMAX)), get_par(par_HRUtype::B_EVAP)))));
//      // if(overFl1 > (prev_Soil+get_dta(tstRM, ts_type::PREF) -evap)) std::cout << "\n problem1 " << get_dta(tstRM, ts_type::PREF) << " d "<< overFl2 << " " << prev_Soil;
//      // if(evap<0) std::cout << "evap " << evap << " prev_Soil - overFl1 " << prev_Soil - overFl1 << " ee " << pow(((get_par(par_HRUtype::SMAX) - next_soil) / get_par(par_HRUtype::SMAX)), get_par(par_HRUtype::B_EVAP)) << std::endl;
//     } else {
//      overFl1 = 0.0;
//      c_init = get_par(par_HRUtype::C_MAX) * (1 - pow((1 - prev_Soil / get_par(par_HRUtype::SMAX)),(1/(get_par(par_HRUtype::B_SOIL) + 1))));
//      next_soil = get_par(par_HRUtype::SMAX) * (1 - pow(1 - (c_init + get_dta(tstRM, ts_type::PREF)) / get_par(par_HRUtype::C_MAX),(get_par(par_HRUtype::B_SOIL) + 1)));
//      overFl2 = std::max(get_dta(tstRM, ts_type::PREF) - get_par(par_HRUtype::SMAX) + next_soil + prev_Soil,0.0);
//      evap =  std::min(static_cast<numberSel>(prev_Soil - overFl2 + get_dta(tstRM, ts_type::PREF) + next_soil), static_cast<numberSel> (get_dta(tstRM, ts_type::PET) * (1 - pow(((get_par(par_HRUtype::SMAX) - next_soil) / get_par(par_HRUtype::SMAX)), get_par(par_HRUtype::B_EVAP)))));
//      // if(overFl2 > (prev_Soil+get_dta(tstRM, ts_type::PREF) -evap)) std::cout << "\n problem " << get_dta(tstRM, ts_type::PREF) << " d "<< overFl2 << " " << prev_Soil << " e " << evap;
//      // if(overFl2>0) std::cout <<"\n" << evap;
//      // if(evap<0) std::cout << "evap " << evap << " prevsoil " << prev_Soil << " ovf "<< overFl2 << " cinit "<< c_init <<" pref "<< get_dta(tstRM, ts_type::PREF) << " next soil "<< next_soil<< " ns-of2 " << (prev_Soil - overFl2 + get_dta(tstRM, ts_type::PREF) + next_soil) << " ee " << pow(((get_par(par_HRUtype::SMAX) - next_soil) / get_par(par_HRUtype::SMAX)), get_par(par_HRUtype::B_EVAP)) << std::endl;
//      // std::cout <<"\nupod1 " << overFl2;
//     }
//
//     // if(prev_Soil<0) {
//       // std::cout <<prev_Soil<<  " 1p_S "<< overFl1 << " ov1 " << overFl2 << " ov2 " << get_dta(tstRM, ts_type::PREF) << " pef "<<  evap << " evap" << "\n";
//     // }
//     if((prev_Soil + get_dta(tstRM, ts_type::PREF))< (evap + overFl1+overFl2)){
//       std::cout  << get_par(par_HRUtype::SMAX)<< " smax "<<prev_Soil<<  " 1p_S "<< overFl1 << " ov1 " << overFl2 << " ov2 " << get_dta(tstRM, ts_type::PREF) << " pef "<<  evap << " evap" << "\n";
//     }
//   prev_Soil = prev_Soil - (overFl1+ overFl2) + get_dta(tstRM, ts_type::PREF) - evap;
//   // if(prev_Soil<0) {
//     // std::cout << get_par(par_HRUtype::SMAX)<< " smax "<< prev_Soil<<  " p_S "<< overFl1 << " ov1 " << overFl2 << " ov2 " << get_dta(tstRM, ts_type::PREF) << " pef "<<  evap << " evap" << "\n";
//   // }
//   set_varValue(prev_Soil, tstRM, ts_type::SOIS);
//   // if(prev_Soil > get_par(par_HRUtype::SMAX)) std::cout << prev_Soil << " efr " << get_dta(tstRM, ts_type::PREF);
//
//   set_varValue(evap,tstRM, ts_type::EVBS);
//   overFL = (overFl1+ overFl2);
//   set_varValue(overFL, tstRM, ts_type::PERC);
//
// // to test
  // aet = get_dta(tstRM,ts_type::EVAC) +  get_dta(tstRM,ts_type::EVAS)  +  evap;
  // set_varValue(aet, tstRM, ts_type::AET);
  // prev_Soil = next_soil;

//
  // Hymod formulation
  c_init = get_par(par_HRUtype::C_MAX) * (1 - pow((1 - prev_Soil / get_par(par_HRUtype::SMAX)),(1/(get_par(par_HRUtype::B_SOIL) + 1))));
  //Overflow if soil tank fully filled
  overFl1 = std::max(static_cast<numberSel>(c_init + get_dta(tstRM, ts_type::PREF) - get_par(par_HRUtype::C_MAX)), static_cast<numberSel>(0.0));
  // if(tstRM ==0){
  //   std::cout << overFl1 <<" prevS "<< prev_Soil << std::endl;
  // }
  ppInf = get_dta(tstRM, ts_type::PREF) - overFl1;
  //Newly proposed soil water depth C
  c_prop = std::min(ppInf + c_init, get_par(par_HRUtype::C_MAX));
    //  //remaining soil input
  //  pref = get_dta(tstRM, ts_type::PREF) -   overFl1;
  //New proposal of state of soil buffer  not affected by evapotranspiration
  next_soil = get_par(par_HRUtype::SMAX) * (1 - pow(1 - c_prop / get_par(par_HRUtype::C_MAX),(get_par(par_HRUtype::B_SOIL) + 1)));
  //Overflow for small C according to Jherman
  overFl2 = std::max(static_cast<numberSel>(0.0),(ppInf - next_soil + prev_Soil));
  //Overflow for small C according to Montanari
  //    overFl2 = std::max(0.0, (c_prop - c_init) - (next_soil - prev_Soil));
  //Evapotranspiration from Soil
  evap =  std::min(static_cast<numberSel>(next_soil), static_cast<numberSel>(get_dta(tstRM, ts_type::PET)*(1 - pow(((get_par(par_HRUtype::SMAX) - next_soil) / get_par(par_HRUtype::SMAX)), get_par(par_HRUtype::B_EVAP)))));
  //Soil buffer state
  next_soil = std::max(static_cast<numberSel>(next_soil - evap),static_cast<numberSel>(0.0));
  //Total overflow
  overFL = overFl1 + overFl2;

set_varValue(next_soil, tstRM, ts_type::SOIS);
set_varValue(evap,tstRM, ts_type::EVBS);
set_varValue(overFL, tstRM, ts_type::PERC);

// aet = get_dta(tstRM,ts_type::EVAC) +  get_dta(tstRM,ts_type::EVAS)  +  evap;
// set_varValue(aet, tstRM, ts_type::AET);

prev_Soil = next_soil;

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
  numberSel prev_Grou1 = 0.0;
  numberSel prev_Grou2 = 0.0;

  switch(_gs_STORtype) {

  case gs_STORtype::LIN_RES:
    BaseOut = prev_Grou * get_par(par_HRUtype::KS) ;
    prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::LINL_RES:
    BaseOut = prev_Grou * get_par(par_HRUtype::KS) ;
    prev_Grou = (prev_Grou * get_par(par_HRUtype::L)) + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::LINBY_RES:
    BaseOut = prev_Grou * get_par(par_HRUtype::KS) + get_par(par_HRUtype::D_BYPASS) * ((1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC));
    prev_Grou = prev_Grou  + (1 - get_par(par_HRUtype::D_BYPASS)) * (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::POW_RES:
    BaseOut = std::pow(prev_Grou, get_par(par_HRUtype::B_EXP)) * get_par(par_HRUtype::KS);
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
    }
    break;

  case gs_STORtype::LIN_2SE:
    BaseOut = prev_Grou2 * get_par(par_HRUtype::KS2);
    BaseOut_1 = prev_Grou1 * get_par(par_HRUtype::KS);

    prev_Grou1 = prev_Grou1 + ((1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC)) - BaseOut_1;
    prev_Grou2 = prev_Grou2 + BaseOut_1 - BaseOut;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou1, tstRM,ts_type::GROS);
    set_varValue(prev_Grou2, tstRM,ts_type::GROS);
    set_varValue(prev_Grou1 + prev_Grou2, tstRM,ts_type::GROS);
    break;

  case gs_STORtype::LIN_2PA:
    BaseOut_1 = prev_Grou1 * get_par(par_HRUtype::KS);
    prev_Grou1 = prev_Grou1 + (1 - get_par(par_HRUtype::ADIV)) * get_par(par_HRUtype::ALPHA) * get_dta(tstRM, ts_type::PERC) - BaseOut_1;

    BaseOut_2 = prev_Grou2 * get_par(par_HRUtype::KS2);
    prev_Grou2 = prev_Grou2 + (1 - get_par(par_HRUtype::ADIV)) * (1 - get_par(par_HRUtype::ALPHA)) * get_dta(tstRM, ts_type::PERC) - BaseOut_2;

    BaseOut = BaseOut_1 + BaseOut_2;

    set_varValue(BaseOut, tstRM, ts_type::BASF);
    set_varValue(prev_Grou1, tstRM,ts_type::GROS);
    set_varValue(prev_Grou2, tstRM,ts_type::GROS);
    set_varValue(prev_Grou1 + prev_Grou2, tstRM,ts_type::GROS);

    break;

  case gs_STORtype::FLEX_RES:
    BaseOut_1 = prev_Grou * get_par(par_HRUtype::KS);
    prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut_1;

    if(get_par(par_HRUtype::THR) > prev_Grou) {
      //lower outlet working
      set_varValue(BaseOut_1, tstRM, ts_type::BASF);

    } else {
      //lower and upper outlets working
      BaseOut_2 = get_par(par_HRUtype::KS2) * (prev_Grou - get_par(par_HRUtype::THR));
      prev_Grou = prev_Grou + (1 - get_par(par_HRUtype::ADIV) ) * get_dta(tstRM, ts_type::PERC) - BaseOut_2;

      BaseOut = BaseOut_1 + BaseOut_2;

      set_varValue(BaseOut, tstRM, ts_type::BASF);
    }

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

  for(ifrb=0; ifrb<help_nmbFR; ifrb++) {
    help_State = get_stateFastres(ifrb);
    helpFastOut = get_par(par_HRUtype::KF) * help_State;
    // help_State = help_State - helpFastOut;
    if(ifrb == 0) {
      help_State = help_State + get_par(par_HRUtype::ADIV) * get_dta(tstRM,ts_type::PERC)- helpFastOut;
    } else help_State = help_State + get_outFastRes((ifrb-1))- helpFastOut;
    set_stateFastRes(help_State,ifrb);
    set_outFastRes(helpFastOut,ifrb);
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
void single_HMunit::interception_NoSnow() {

  numberSel CanOut = 0.0, StemOut = 0.0, OverflowCan = 0.0, OverflowStem, EvapCanop = 0.0, EvapStem = 0.0, Througf = 0.0;

  OverflowCan = std::max((prevCanS - get_par(par_HRUtype::CAN_ST)),0.0);
  //!< VIC model for canopy evaporation (prevCanS/ get_par(par_HRUtype::CAN_ST))^(2/3)
  prevCanS = prevCanS - OverflowCan;
  EvapCanop = std::min(std::pow(((prevCanS) / get_par(par_HRUtype::CAN_ST)),2/3),prevCanS);
  prevCanS = prevCanS - EvapCanop;

  CanOut = std::min((prevCanS / get_par(par_HRUtype::CAN_ST) * EvapCanop),prevCanS);
  prevCanS = prevCanS - CanOut;
  set_varValue(prevCanS, tstRM, ts_type::CANS);

  prevCanS =  prevCanS + get_par(par_HRUtype::CDIV) * (get_dta(tstRM, ts_type::PREC) + get_dta(tstRM, ts_type::MELT));

  //!< VIC model for canopy evaporation (prevCanS/ get_par(par_HRUtype::CAN_ST))^(2/3)
  EvapStem = std::min(std::pow(((prevSteS) / get_par(par_HRUtype::STEM_ST)),(2/3)), prevSteS);
  prevSteS = prevSteS - EvapStem;

  OverflowStem = std::max((prevSteS - get_par(par_HRUtype::STEM_ST)),0.0);
  prevSteS = prevSteS - OverflowStem;

  StemOut = std::min((prevCanS) / get_par(par_HRUtype::CAN_ST) * EvapStem, prevSteS);
  prevSteS = prevSteS - StemOut;
  set_varValue(prevSteS, tstRM, ts_type::STES);

  prevSteS = prevSteS + get_par(par_HRUtype::SDIV) * (get_dta(tstRM, ts_type::PREC) + get_dta(tstRM, ts_type::MELT)) + (1 - get_par(par_HRUtype::CSDIV)) * (CanOut + OverflowCan);

  set_varValue((CanOut + OverflowCan), tstRM, ts_type::CANF);
  set_varValue(EvapCanop, tstRM, ts_type::EVAC);

  set_varValue((StemOut + OverflowStem), tstRM, ts_type::STEF);
  set_varValue(EvapStem, tstRM, ts_type::EVAS);

  Througf = (OverflowCan + CanOut) * get_par(par_HRUtype::CSDIV) + StemOut + OverflowStem;

  set_varValue(Througf, tstRM, ts_type::TROF);

  set_varValue((get_dta(tstRM, ts_type::CANS) + get_dta(tstRM, ts_type::STES)),tstRM,ts_type::INTS);

  return ;
}



/** \brief updates states in interception storage in single HRU unit
 *
 *  updating of states in canopy and stem-trunk storage in the single hru unit
 *  calculating related fluxes
 *
 */
void single_HMunit::interception_WithSnow() {
  //  numberSel CanOut = 0.0, StemOut = 0.0, OverflowCan = 0.0, OverflowStem, EvapCanop = 0.0, EvapStem = 0.0, Througf = 0.0;
  numberSel OverflowCan = 0.0, OverflowStem= 0.0, EvapCanop = 0.0, EvapStem = 0.0, Througf = 0.0;

  OverflowCan = std::max((prevCanS - get_par(par_HRUtype::CAN_ST)),0.0);
  prevCanS = prevCanS - OverflowCan;
  //!< VIC model for canopy evaporation (prevCanS/ get_par(par_HRUtype::CAN_ST))^(2/3)
   EvapCanop = std::min(pow(((prevCanS) / get_par(par_HRUtype::CAN_ST)),2/3),prevCanS);
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
    interception_WithSnow();
    //    std::cout << " snow " << Snoww << " \n";
  } else {
    Snoww = 0.0;
    set_varValue(Snoww,tstRM,ts_type::SNOW);
    snow_melt();
    interception_NoSnow();
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
  if( !(Numdta >0)) {
    std::cout << std::endl << "There is an error in data loadings." << std::endl;
    std::cout << "It is impossible to calculate the basic for loop in runHB function" << std:: endl << "It is controlled by the length " << Numdta << "." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  numberSel helprm=0.0;
  for(tstRM=0; tstRM < Numdta ; tstRM++) {
    interception_snow();
    surface_retention();
    soil_buffer();
    slow_response(gs_STORAGE);
    fast_response();
    helprm = (get_dta(tstRM,ts_type::BASF) + get_dta(tstRM,ts_type::DIRR));
    set_varValue(helprm ,tstRM,ts_type::TOTR);

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

  if(helpnumDTA != temp_input.size()) {
    std::cout << "Different number of time intervals in precipitation input " << prec_input.size() \
              << " the temperature input has " << temp_input.size() << " inputs." << std::endl;
    std::exit(EXIT_FAILURE);

  }
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

  } else {

    std::cout << "There is a error in opening and reading from the file with path and name: " << Filet.c_str() << std::endl;

  }

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

/** \brief Prints all names and values of parameters of HRU
 *
 *
 */
void single_HMunit::print_Pars() {

  std::cout << "\nThe HRU ID " << getIdHru() << std::endl;

  par_HRU.p_param();

  return ;

}

void single_HMunit::set_GStype(gs_STORtype _gs_STORtype) {

  gs_STORAGE = _gs_STORtype;
  print_GStype();

  return ;

}

gs_STORtype single_HMunit::get_GStype() {
  return gs_STORAGE;

}

void single_HMunit::print_GStype() {

  switch(gs_STORAGE) {

  case gs_STORtype::LIN_RES:

    std::cout << "The gs_STORE is a LIN reservoir." << std::endl;

    break;

  case gs_STORtype::LINL_RES:

    std::cout << "The gs_STORE is a LINL reservoir." << std::endl;

    break;

  case gs_STORtype::LINBY_RES:

    std::cout << "The gs_STORE is a LINBY reservoir." << std::endl;

    break;

  case gs_STORtype::LIN_2SE:

    std::cout << "The gs_STORE is a LIN_2SE reservoir." << std::endl;

    break;

  case gs_STORtype::LIN_2PA:

    std::cout << "The gs_STORE is a LIN_2PA reservoir." << std::endl;

    break;

  case gs_STORtype::POW_RES:

    std::cout << "The gs_STORE is a POW reservoir." << std::endl;

    break;

  case gs_STORtype::EXP_RES:

    std::cout << "The gs_STORE is a EXP reservoir." << std::endl;

    break;

  case gs_STORtype::FLEX_RES:

    std::cout << "The gs_STORE is a FLEX reservoir." << std::endl;

    break;

  }

}

