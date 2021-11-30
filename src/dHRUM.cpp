#include "dHRUM.h"

dHRUM::dHRUM(): dHruVec(),
  dHruVecId(),
  dimHM(0),
  basinArea(0),
  Areas(0),
  numTs(0),
  numParsAllHRus(0),
  basinDta(),
  gs_STORtypes(),
  NumFastRes(0){
  //ctor
}

dHRUM::~dHRUM() {
  //dtor
}

dHRUM::dHRUM(const dHRUM& other): dHruVec(),
  dHruVecId(),
  dimHM(0),
  basinArea(0),
  Areas(0),
  numTs(0),
  numParsAllHRus(0),
  basinDta(),
  gs_STORtypes(),
  NumFastRes(0){

  dHruVec = other.dHruVec;
  dHruVecId = other.dHruVecId;
  dimHM = other.dimHM;
  basinArea = other.basinArea;
  Areas = other.Areas;
  numTs = other.numTs;
  numParsAllHRus = other.numParsAllHRus;
  basinDta = other.basinDta;
  gs_STORtypes = other.gs_STORtypes;
  NumFastRes = other.NumFastRes;
}

dHRUM& dHRUM::operator=(const dHRUM& rhs) {
  if (this == &rhs) return *this;
  else {
    dHruVec = rhs.dHruVec;
    dHruVecId = rhs.dHruVecId;
    dimHM = rhs.dimHM;
    basinArea = rhs.basinArea;
    Areas = rhs.Areas;
    numTs = rhs.numTs;
    numParsAllHRus = rhs.numParsAllHRus;
    basinDta = rhs.basinDta;
    gs_STORtypes = rhs.gs_STORtypes;
    NumFastRes = rhs.NumFastRes;
  }
  return *this;
}

/** \brief Initialization of vectors of single HM units
 *
 * \param number of single HM units units
 * \param instance of single HM units unit
 *
 */
void dHRUM::initHrusVec(const unsigned& numHMunits, const single_HMunit& shru) {

  for(unsigned int i=0; i<numHMunits; i++)  {
    dHruVec.push_back(shru);
  }

  dimHM = numHMunits;

  return ;

}


void dHRUM::setInputsToAllHrus(std::string namesFilePath) {
  //  std::cout << "\n dssda " << dimHM ;

  //  #pragma omp parallel
  //  {
#pragma omp parallel for
  for(unsigned it=0; it<dimHM; it++) {
                    // std::cout << "Loading the input data  to HRU ID " << it << std::endl;
    //    dHruVec[it].set_paramsToSim(parsToLoad);

    dHruVec[it].read_InputFromFile(namesFilePath.c_str());
    //  dHruVec[it].set_PetVars(50.1,pet_Type::OUDIN);
    //  dHruVec[it].calc_Pet();
    //  dHruVec[it].run_HB();
    // std::cout << "Loading the input data  to HRU ID " << it << std::endl;
  }
  //  }
  initdHRUbasinDTA();
  //  std::cout << " Number of TS " << numTs << std::endl;
  return ;
}

void dHRUM::setInputsToOneHru(std::string namesFilePath, unsigned Id) {

  dHruVec[Id].read_InputFromFile(namesFilePath.c_str());

}

void dHRUM::initGWtypeToAlldHrus(std::vector<std::pair<unsigned,gs_STORtype>>& gs_STORtypes) {

  for(unsigned int i=0; i<gs_STORtypes.size(); i++) {
    dHruVec[gs_STORtypes[i].first].set_GStype(gs_STORtypes[i].second);
  }

  return;

}


void dHRUM::initSoilStypeToAlldHrus(std::vector<std::pair<unsigned,soil_STORtype>>& soil_STORtypes) {

    for(unsigned int i=0; i<soil_STORtypes.size(); i++) {
      dHruVec[soil_STORtypes[i].first].set_soilStorType(soil_STORtypes[i].second);
    }

    return;

  }

void dHRUM::setParamsToAlldHrus(std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad) {
  //  #pragma omp parallel
  //  {
#pragma omp parallel for
  for(unsigned it=0; it<dimHM; it++) {
    dHruVec[it].set_paramsToSim(parsToLoad);
    //    dHruVec[it].read_InputFromFile(namesFilePath.c_str());
    //  dHruVec[it].set_PetVars(50.1,pet_Type::OUDIN);
    //  dHruVec[it].calc_Pet();
    //  dHruVec[it].run_HB();
    //    std::cout << "Setting the params to the input data to HRU ID " << it << std::endl;
  }
  //  }
  return ;
}

std::vector<std::string> dHRUM::getRequiredParamsForHru(unsigned hruId) {

  gs_STORtype gwType = dHruVec[hruId].get_GStype();
  //soil_STORtype swType = dHruVec[hruId].get_soilStorType();

  std::vector<std::string> params = {};
  //std::vector<std::string> params2 = {};
  // params.reserve(50);

  switch(gwType) {
    case gs_STORtype::LIN_RES:
      params = {"KS","ADIV"};
      break;
    case gs_STORtype::LINL_RES:
      params = {"KS","L","ADIV"};
      break;
    case gs_STORtype::LINBY_RES:
      params = {"KS","D_BYPASS","ADIV"};
      break;
    case gs_STORtype::POW_RES:
      params = {"KS","B_EXP","ADIV"};
      break;
    case gs_STORtype::EXP_RES:
      params = {"KS","B_EXP","ADIV"};
      break;
    case gs_STORtype::LIN_2SE:
      params = {"KS","KS2","ADIV"};
      break;
    case gs_STORtype::LIN_2PA:
      params = {"KS","KS2","ALPHA","ADIV"};
      break;
    case gs_STORtype::FLEX_RES:
      params = {"KS","KS2","THR","ADIV"};
      break;
  }

  // switch(swType) {
  // case soil_STORtype::PDM:
  //   params2 = {"SMAX","CMIN","C_MAX","B_SOIL","B_EVAP"};
  //   break;
  // case soil_STORtype::COLLIE_V2:
  //   params2 = {"SMAX","FOREST_FRACT","FC","KF"};
  //   break;
  // case soil_STORtype::NEW_ZEALAND:
  //   params2 = {"FC","FOREST_FRACT","FC","SMAX","KF","KF_NONLIN","KF2"};
  //   break;
  // case soil_STORtype::GR4J:
  //   params2 = {"SMAX"};
  //   break;
  // case soil_STORtype::SBROOK_V1:
  //   params2 = {"SMAX","FOREST_FRACT","FC","KF_NONLIN"};
  //   break;
  // case soil_STORtype::HILLSLOPE:
  //   params2 = {"KF_NONLIN","SMAX","C"};
  //   break;
  // case soil_STORtype::PLATEAU:
  //   params2 = {"INFR_MAX","WP","SMAX","RF","C"};
  //   break;
  // }

  //std::vector<std::string> v(params);
  //
  //v.insert(v.end(), params2.begin(), params2.end());

  // for(unsigned i=0;i<v.size();i++){
  //   std::cout << v[i] <<  std::endl;
  // }

  return params;
}

void dHRUM::setParamsToOneHru(std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad, unsigned Id) {

  dHruVec[Id].set_paramsToSim(parsToLoad);

  return ;

}

void dHRUM::calcPetToAllHrus(numberSel Latit, pet_Type myPetType) {
  //
  //  std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad;
  //
  //  parsToLoad.push_back(std::make_pair(1.23, par_HRUtype::B_SOIL));
  //  parsToLoad.push_back(std::make_pair(500.0, par_HRUtype::C_MAX));
  //  parsToLoad.push_back(std::make_pair(1.30, par_HRUtype::B_EVAP));
  //  parsToLoad.push_back(std::make_pair(0.20, par_HRUtype::KS));
  //  parsToLoad.push_back(std::make_pair(0.40, par_HRUtype::KF));
  //  parsToLoad.push_back(std::make_pair(0.60, par_HRUtype::ADIV));
  //  parsToLoad.push_back(std::make_pair(0.60, par_HRUtype::CDIV));
  //  parsToLoad.push_back(std::make_pair(0.30, par_HRUtype::SDIV));
  //  parsToLoad.push_back(std::make_pair(2.00, par_HRUtype::CAN_ST));
  //  parsToLoad.push_back(std::make_pair(2.00, par_HRUtype::STEM_ST));
  //  parsToLoad.push_back(std::make_pair(0.50, par_HRUtype::CSDIV));
  //  parsToLoad.push_back(std::make_pair(0.50, par_HRUtype::DDFA));
  //  parsToLoad.push_back(std::make_pair(2.00, par_HRUtype::RETCAP));
  //

  //  #pragma omp parallel
  //  {
#pragma omp parallel for
  for(unsigned it=0; it<dimHM; it++) {
    ////    dHruVec[it].set_paramsToSim(parsToLoad);
    //    dHruVec[it].read_InputFromFile(namesFilePath.c_str());
    dHruVec[it].set_PetVars(Latit,myPetType);
    dHruVec[it].calc_Pet();
    //  dHruVec[it].run_HB();
    //    std::cout <<  "Calculating the PET data to HRU ID " << it << std::endl;
  }
  //  }
  return ;

}

void dHRUM::calcPetToOneHru(numberSel Latit, pet_Type myPetType, unsigned HruId){

  dHruVec[HruId].set_PetVars(Latit,myPetType);
  dHruVec[HruId].calc_Pet();

  return ;

}


void dHRUM::calcPetToAllHrusDist(hdata LatitVec, std::vector<std::string> petTypeString){

  std::map<std::string, pet_Type> s_mapString_ToPet_Type = {
    {"OUDIN", pet_Type::OUDIN},
    {"HAMON", pet_Type::HAMON },
    {"THORNTHWAITE",  pet_Type::THORNTHWAITE},
    {"BLANEYCRIDDLE",  pet_Type::BLANEYCRIDDLE},
    {"JENSENHAISE",   pet_Type::JENSENHAISE},
    {"MCGUINNESSBORDNE",  pet_Type::MCGUINNESSBORDNE}
  };

#pragma omp parallel for
  for(unsigned it=0; it<dimHM; it++) {
    switch(s_mapString_ToPet_Type[petTypeString[it]]) {
     case pet_Type::OUDIN:
      calcPetToOneHru(LatitVec[it], pet_Type::OUDIN, it);
      break;
    case pet_Type::HAMON:
      calcPetToOneHru(LatitVec[it], pet_Type::HAMON, it);
      break;
    case pet_Type::THORNTHWAITE:
      calcPetToOneHru(LatitVec[it], pet_Type::THORNTHWAITE, it);
      break;
    case pet_Type::BLANEYCRIDDLE:
      calcPetToOneHru(LatitVec[it], pet_Type::BLANEYCRIDDLE, it);
      break;
    case pet_Type::JENSENHAISE:
      calcPetToOneHru(LatitVec[it], pet_Type::JENSENHAISE, it);
      break;
    case pet_Type::MCGUINNESSBORDNE:
      calcPetToOneHru(LatitVec[it], pet_Type::MCGUINNESSBORDNE, it);
      break;
         }
  }

  return ;

}

void dHRUM::calcHbToAllHrus() {
  //
  //  std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad;
  //
  //  parsToLoad.push_back(std::make_pair(1.23, par_HRUtype::B_SOIL));
  //  parsToLoad.push_back(std::make_pair(500.0, par_HRUtype::C_MAX));
  //  parsToLoad.push_back(std::make_pair(1.30, par_HRUtype::B_EVAP));
  //  parsToLoad.push_back(std::make_pair(0.20, par_HRUtype::KS));
  //  parsToLoad.push_back(std::make_pair(0.40, par_HRUtype::KF));
  //  parsToLoad.push_back(std::make_pair(0.60, par_HRUtype::ADIV));
  //  parsToLoad.push_back(std::make_pair(0.60, par_HRUtype::CDIV));
  //  parsToLoad.push_back(std::make_pair(0.30, par_HRUtype::SDIV));
  //  parsToLoad.push_back(std::make_pair(2.00, par_HRUtype::CAN_ST));
  //  parsToLoad.push_back(std::make_pair(2.00, par_HRUtype::STEM_ST));
  //  parsToLoad.push_back(std::make_pair(0.50, par_HRUtype::CSDIV));
  //  parsToLoad.push_back(std::make_pair(0.50, par_HRUtype::DDFA));
  //  parsToLoad.push_back(std::make_pair(2.00, par_HRUtype::RETCAP));
  //
  //  #pragma omp parallel
  //  {
#pragma omp parallel for

  for(unsigned it=0; it<dimHM; it++) {
    ////    dHruVec[it].set_paramsToSim(parsToLoad);
    //    dHruVec[it].read_InputFromFile(namesFilePath.c_str());
    //  dHruVec[it].set_PetVars(50.1,pet_Type::OUDIN);
    //  dHruVec[it].calc_Pet();
    dHruVec[it].run_HB();
    //    std::cout << "Calculating the HB data to HRU ID " << it << " ";
  }
  //  }
  return ;

}

void dHRUM::setAreasToHrus(hdata vec_Areas) {

  Areas = vec_Areas;
  //       std::cout << "Areas " << Areas[0] << std::endl;
  //
  //  #pragma omp parallel
  //  {
#pragma omp parallel for

  for(unsigned it=0; it<dimHM; it++) {
    dHruVec[it].set_Area(Areas[it]);
    //     std::cout << "HRu Area " << dHruVec[it].get_Area() << std::endl;
  }
  //  }
  setBasinArea();

  return ;
}


void dHRUM::setBasinArea() {

  basinArea = Areas.sum();

  return ;

}

void dHRUM::gatherTsFromHrus() {
  //ToDo parallel challenge
  for(const auto& it : all_ts) {
    basinDta.setOneTstoZero(it);
  }
  //  for(const auto& itr : all_ts) {
#pragma omp parallel for
  for(unsigned int itRts=0; itRts<numTSvars; itRts++) {
    ts_type    itTS =all_ts[itRts];
    hdata helpValAr(0.0,numTs);
    //    std::cout << helpValAr[0] << " " << dHruVec[0].getSingleHruTsDta(it)[0] << " fdsoi " <<std::endl;
    //  dHruVec[itHru].getSingleHruTsDta(it) * Areas[itHru] /
    // dHruVec[0].getSingleHruTsDta(it) *
    for(unsigned itHru=0; itHru<dimHM; itHru++) {
      helpValAr +=   dHruVec[itHru].getSingleHruTsDta(itTS) * Areas[itHru] / basinArea;
      //      std::cout << helpValAr[0] << " " << Areas[0] << " ba " << basinArea <<std::endl;
      //      for(unsigned ii=0;ii<numTs;ii++) std::cout << helpValAr[ii] << "\n";
    }
    basinDta.s_data(helpValAr,itTS,false) ;
    //    std::cout << helpValAr[0] << " uiofsdu " << dHruVec[0].getSingleHruTsDta(itTS)[0] <<std::endl;
    //    dHruVec[0].getSingleHruTsDta(it);
  }

  return ;

}


numberDta dHRUM::get_numTS() {

  return numTs;

}


void dHRUM::set_numTS() {

  numTs = dHruVec[0].get_numdta();

  return ;

}

void dHRUM::initdHRUbasinDTA() {

  set_numTS();
  basinDta = dHruVec[0].getAllData();
  //  printAllDta("./data/tests/outALL/outALLdta.txt");
  basinDta.setAllToZeros(false,true);

  return ;

}

void dHRUM::printAllDta(const std::string& Filet) {

  basinDta.printDataToFile(Filet);

  return ;

}


void dHRUM::initHrusID(const std::vector<std::string>& SingleHrusIds) {

  for(unsigned ii=0;ii<dimHM;ii++){
    dHruVec[ii].setIdHru(SingleHrusIds[ii]);
  }

  return ;

}

std::string dHRUM::getSingleHruId(unsigned hruId) {

  return dHruVec[hruId].getIdHru();

}

hdata dHRUM::get_HbDta(const ts_type& _tsType){

  return basinDta.get_HbTsData(_tsType);

}


caldata dHRUM::get_CalDta(const cal_Type& _calType){

  return basinDta.getCalData(_calType);

}


void dHRUM::loadPTDatToAllHrus(hdata Prec, hdata Temp, const numberSel& val,const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay) {

#pragma omp parallel for
  for(unsigned it=0; it<dimHM; it++) {
    dHruVec[it].load_data_PT(Prec,Temp,val,inYear,inMonth,inDay);
  }

  initdHRUbasinDTA();

  return ;
}

numberDta dHRUM::getdHRUdim() {

  return dimHM;

}

void dHRUM::load_CalDataToAllHrus(const caldata& yyear, const caldata& mmonth, const caldata& dday) {

#pragma omp parallel for
 for(unsigned it=0; it<dimHM; it++) {
   dHruVec[it].load_calData(yyear, mmonth, dday);
 }

 return ;

}


void dHRUM::load_PrecTempToAllHrus(const hdata& Prec, const hdata& Temp) {

  unsigned numDTA = 0;
  numDTA = Prec.size();
// std::cout << numDTA << "\n";

#pragma omp parallel for
  for(unsigned it=0; it<dimHM; it++) {
    dHruVec[it].init_inputs(0.0,numDTA);
    dHruVec[it].set_data_prec_temp(Prec, Temp);
  }

  return ;

}
unsigned dHRUM::get_singleHRUnumPars(unsigned Id){

  return (dHruVec[Id].get_numPars());

}

void dHRUM::set_numPars() {

  numParsAllHRus.resize(dimHM);

#pragma omp parallel for
  for(unsigned it=0;it<dimHM;it++){
    numParsAllHRus[it] = get_singleHRUnumPars(it);
    //std::cout <<"The number of params for HRU unit "<< it<< " is " << get_singleHRUnumPars(it) << std::endl;
  }

  return;

}

void dHRUM::print_Pars() {

  for(unsigned it; it <dimHM;it++){
    dHruVec[it].print_Pars();
  }

  return;

}


void dHRUM::loadPTInputsToOneHru(hdata Prec, hdata Temp, const numberSel& val,const unsigned& inYear, const unsigned& inMonth,const unsigned& inDay, unsigned HruIt) {

  dHruVec[HruIt].load_data_PT(Prec,Temp,val,inYear,inMonth,inDay);

  return;

}

std::vector<std::string> dHRUM::getHRUIds() {

  std::vector<std::string> ids(getdHRUdim());
  for(unsigned int i=0; i<ids.size(); i++) {
    ids[i] = getSingleHruId(i);
  }

  return ids;
}

numberSel dHRUM::getTsDta(const ts_type& _tsType, const unsigned& HruIndex, const unsigned& tst){

  return dHruVec[HruIndex].get_dta(tst, _tsType);

}


void dHRUM::set_numFastReservoirsToHrus(){

#pragma omp parallel for
  for(unsigned it=0;it<dimHM;it++){
    dHruVec[it].set_nmbFastres(NumFastRes[it]);
  }

}


void dHRUM::set_numFastReservoirs(caldata numFR){
  unsigned numNFR=0;

  NumFastRes.resize(dimHM);

#pragma omp parallel for
  for(unsigned it=0;it<dimHM;it++){
    NumFastRes[it] = numFR[it];
  }

}
