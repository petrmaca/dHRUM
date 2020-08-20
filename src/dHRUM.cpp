#include "dHRUM.h"

dHRUM::dHRUM(): dHruVec(),
  dimHM(0),
  basinArea(0),
  Areas(0),
  numTs(0),
  numParsAllHRus(0),
  basinDta() {
  //ctor
}

dHRUM::~dHRUM() {
  //dtor
}

dHRUM::dHRUM(const dHRUM& other): dHruVec(),
  dimHM(0),
  basinArea(0),
  Areas(0),
  numTs(0),
  numParsAllHRus(0),
  basinDta() {

  dHruVec = other.dHruVec;
  dimHM = other.dimHM;
  basinArea = other.basinArea;
  Areas = other.Areas;
  numTs = other.numTs;
  numParsAllHRus = other.numParsAllHRus;
  basinDta = other.basinDta;
}

dHRUM& dHRUM::operator=(const dHRUM& rhs) {
  if (this == &rhs) return *this;
  else {
    dHruVec = rhs.dHruVec;
    dimHM = rhs.dimHM;
    basinArea = rhs.basinArea;
    Areas = rhs.Areas;
    numTs = rhs.numTs;
    numParsAllHRus = rhs.numParsAllHRus;
    basinDta = rhs.basinDta;
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


void dHRUM::setParamsToAllHrus(std::vector<std::pair<numberSel,par_HRUtype>> parsToLoad) {
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
    std::cout <<"The number of params for HRU unit "<< it<< " is " << get_singleHRUnumPars(it) << std::endl;
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
