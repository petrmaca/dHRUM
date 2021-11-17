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
  void interception_NoSnow();//!< Update the Canopy and Stem Interception storages without snow
  void interception_WithSnow();//!< Update the Canopy and Stem Interception storage with snow
  void surface_retention();//!< Update surface retention
  void soil_buffer();//!< Update the soil buffer states
  void fast_response();//!< The fast runoff response
  void slow_response(gs_STORtype _gs_STORtype);//!< The slow runoff response

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

  numberSel et_demand;//!< The helper on ET demand

  numberDta help_nmbFR;//!< The helper for number of fast reservoirs
  numberDta ifrb;//!< For loop counter

  numberSel Area;//!< The area of HM unit in m2
  std::string IdHru;

  gs_STORtype gs_STORAGE;

};

#endif // single_HMunit_H
