#ifndef SINGLE_HMUNIT_H
#define SINGLE_HMUNIT_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include "numberSel.h"
#include "pondSel.h"
#include "data_HB_1d.h"
#include "params.h"


class single_HMunit {
 public:
  single_HMunit();
  ~single_HMunit();
  single_HMunit(const single_HMunit& other);
  single_HMunit& operator=(const single_HMunit& other);

  numberSel get_dta(const unsigned& tst, const ts_type& _tsType);//!< Getter data on HB of single pdm unit
  void set_data_prec_temp(const hdata& _prec_dta,const hdata& _temp_dta);//!< Setter data on HB of single pdm unit
  void set_data(const hdata& dta,const ts_type&_tsType);//!< Setter data on HB of single pdm unit
  void set_varValue(const numberSel& dta,const  unsigned& tst,const ts_type& _tsType);

  void set_calender();
  void calc_Pet();
  void set_PetVars(const numberSel& newLatitude, const pet_Type& newPeType);

  void snow_melt();//!< The update opf snow storage and estimating the snowmelt
  void interception_snow();//!< Update snow storage and interception
  void interception_NoSnow(interception_STORtype _intrc_STORAGE);//!< Update the Canopy and Stem Interception storages without snow
  void interception_WithSnow(interception_STORtype _intrc_STORAGE);//!< Update the Canopy and Stem Interception storage with snow
  void surface_retention(surface_STORtype _surf_STORtype);//!< Update surface retention
  void soil_buffer(soil_STORtype _soil_STORtype);//!< Update the soil buffer states
  void fast_response(fast_Response _fast_RESPONSE);//!< The fast runoff response
  void slow_response(gs_STORtype _gs_STORtype);//!< The slow runoff response
  void upadate_actualET();//!< The estimation of actual evapotranspiration

  void run_HB();//!< Update hydrological balance (HB) of single pdm unit run single pdm model

  void set_params();//!< Setter on params of single pdm unit
  numberSel get_par(const par_HRUtype& _parType);//!< Getter of single one parameter on single pdm unit
  unsigned get_numPars();//!< Getter for number of parameters

  void set_ZeroinitStates(const unsigned& numres);//!< The setter of initial states of state variables
  void set_ZeroStates();//!< The setter of zero states of state variables
  numberSel get_initState(const init_Stype& _Stype);//!< Getter for initial state of given state variable

  unsigned get_numdta();//!< Get the total number of intervals in time series
  void set_nmbFastres(const unsigned& nmbRes);//!< Set number of fast linear reservoirs
  unsigned get_nmbFastRes();
  numberSel get_outFastRes(const unsigned& itFasRes);
  numberSel get_stateFastres(const unsigned& itFasRes);
  void set_outFastRes(const numberSel& helpOut,const unsigned& itFasRes);
  void set_stateFastRes(const numberSel& helpState,const unsigned& itFasRes);

  void init_inputs(numberSel val, unsigned numDTA);
  void load_data_PT(const hdata& prec_input, const hdata& temp_input, const numberSel& val,const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay);
  void set_calender(const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay,const unsigned& initNumTS);
  void load_calData(const caldata& yyear, const caldata& mmonth, const caldata& dday);

  void set_paramsToSim(std::vector<std::pair<numberSel,par_HRUtype>>& parsToLoad);
  void p_defaultParams(bool Bprint);

  void print_OutputToFile(const std::string& Filet);
  void read_InputFromFile(const std::string& Filet);

  void set_Area(numberSel area);
  numberSel get_Area();

  hdata getSingleHruTsDta(const ts_type& _tsType);
  data_HB_1d getAllData();

  void setIdHru(const std::string& IdToSet);
  std::string getIdHru();
  void print_Pars();

  void set_GStype(gs_STORtype _gs_STORtype);
  gs_STORtype get_GStype();
  // void print_GStype();

  numberSel update_ETDEMAND(const numberSel& ET, bool ET_demand);

  void set_soilStorType(soil_STORtype _soil_STORtype);
  soil_STORtype get_soilStorType();
  void print_soilStorType();

  void set_inteceptionType(interception_STORtype _intrc_STORAGE);
  interception_STORtype get_intercetionStorType();
  // void print_soilStorType();

  void set_surfaceStor(surface_STORtype _srfs_STORAGE);
  surface_STORtype get_surfaceStorType();

  void set_fast_response(fast_Response _fast_RESPONSE);
  fast_Response get_fastResponseType();

  std::vector<numberSel> water_balance(numberSel next_soil, numberSel val, std::vector<numberSel> vals);//!< Method for preserving mass balance

  void print_sHRU_settings();
  void print_gs_STORtype();
  void print_interception_STORtype();
  void print_surface_STORtype();
  void print_fastresponseType();
  void print_pondType();


  void current_params();

  std::vector<double>Current_par_val;
  std::vector<double>Current_uppar_val;
  std::vector<double>Current_lowpar_val;

  void ponds(pond_type _pondtype);
  void set_pond_type(pond_type _pondtype);
  pond_type get_pondtype();
  void set_pond_variables(std::vector<std::pair<std::string,numberSel>>& PondDefs,std::vector<std::pair<std::string,std::string>>& PondBeh);

  //pond switches
  numberSel pond_ET(ETpond_type _etpond_type);
  numberSel pond_SOISperc(PondSOISPerc_type _soispond_type);
  numberSel pond_GWperc(PondGWPerc_type _gwpond_type);
  numberSel pond_regular_out(PondRouT_type _RouT_type);


protected:

private:
  numberDta tstRM;//!< The counter for main loop in run model
  params par_HRU;//!< The parameters in PDM instances

  data_HB_1d hyd_dta;//!< The data of all time series of hydrological variables

  numberSel prev_Soil;//!< The helper variable for updating soil storage
  numberSel prev_Grou;//!< The helper variable for updating groundwater storage
  numberSel prevCanS;//!<  The helper variable for Canopy interception storage
  numberSel prevSteS;//!<  The helper variable for Stem interception storage
  numberSel prevSnoS;//!<  The helper variable for Snow storage
  numberSel prev_SurS;//!< The helper variable for updating surface storage

  numberSel prev_GroS1;//!< The helper variable for updating groundwater storage for LIN_2SE and LIN_2PA
  numberSel prev_GroS2;//!< The helper variable for updating groundwater storage for LIN_2SE and LIN_2PA

  numberSel et_demand;//!< The helper on ET demand

  numberDta help_nmbFR;//!< The helper for number of fast reservoirs
  numberDta ifrb;//!< For loop counter

  numberSel Area;//!< The area of HM unit in m2
  std::string IdHru;



  gs_STORtype gs_STORAGE;//!< Type of groundwater storage
  soil_STORtype soil_STORAGE;//!< Type of soil storage
  interception_STORtype intrc_STORAGE;//!< type of iterception storage
  surface_STORtype srfs_STORAGE;//!< the type of surface retentions storage
  fast_Response fast_RESPONSE;//!< the type of fast response
  pond_type pond;//!< the type of pond (water reservoir)
  ETpond_type ET_POND;//!< Type of evaporation from pond surface
  PondSOISPerc_type pondSOISPERCin;//!< Type of percolation from   soil to pond
  PondSOISPerc_type pondSOISPERCout;//!< Type of percolation from pond to soil
  PondGWPerc_type  pondGWPERCin;//!< Type of percolation from   groundwater to pond
  PondGWPerc_type  pondGWPERCout;//!< Type of percolation from pond to groundwater
  PondRouT_type  PondROUT ;//!< Pond outlet method


  numberSel pondArea;  //!< The area of the pond [m2]
  numberSel PonsMax; //!< The maximum pond volume [m3]
  numberSel MRF; //!< Minimum residual flow (MZP) [m3/s]
  numberSel Coflw; //!< Constant user defined regular outflow (RouT) from pond [m3/s]


  std::vector<std::string>Current_par_names;


};

#endif // single_HMunit_H
