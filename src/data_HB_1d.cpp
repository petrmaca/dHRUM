#include "data_HB_1d.h"

data_HB_1d::data_HB_1d(): numTS(0),
  Latitude(0.0),
  PETtype(pet_Type::OUDIN),
  year(1,1),
  month(1,1),
  day(1,1),
  Jday(1,1),
  init_year(1,1),
  init_month(1,1),
  init_day(1,1),
  Prec(1,1),
  Snow(1,1),
  AEt(1,1),
  PEt(1,1),
  Temp(1,1),
  TroF(1,1),
  SteF(1,1),
  CanF(1,1),
  CanS(1,1),
  SteS(1,1),
  EvaC(1,1),
  EvaS(1,1),
  EvbS(1,1),
  IntS(1,1),
  SoiS(1,1),
  GroS(1,1),
  SurS(1,1),
  TotR(1,1),
  Basf(1,1),
  DirR(1,1),
  Melt(1,1),
  Perc(1,1),
  Pref(1,1),
  Etsw(1,1),
  PonS(1,1),
  EtpO(1,1),
  PoiS(1,1),
  PoiG(1,1),
  init_SoiS(0.0),
  init_GroS(0.0),
  init_CanS(0.0),
  init_SteS(0.0),
  init_SnoS(0.0),
  init_SurS(0.0),
  init_GroS1(0.0),
  init_GroS2(0.0),
  numfastRes(1),
  StateFastRes(1.0,1),
  OutFastRes(1.0,1) {
  //ctor
//
//  Prec.resize(1, 999999);
//
//  AEt = Prec;
//  PEt = Prec;
//  Snow = Prec;
//  Temp = Prec;
//
//  TroF = Prec;
//  SteF = Prec;
//  CanF = Prec;
//
//  EvaC = Prec;
//  EvaS = Prec;
//  EvbS = Prec;
//
//  CanS = Prec;
//  SteS = Prec;
//
//  SoiS = Prec;
//  GroS = Prec;
//  SurS = Prec;
//  IntS = Prec;
//  TotR = Prec;
//  Basf = Prec;
//  DirR = Prec;
//  Melt = Prec;
//  Perc = Prec;
//  Pref = Prec;
//
//
//  numTS = Prec.size();
//
//  std::cout << "\nSize of Prec : " << Prec.size() << "\n";
//  std::cout << "\nSize of Prec : " << Prec.size() << " " << Prec[0] << "\n";
//
//  if(Prec.size()!= Perc.size()) {
//    std::cout << "\n Something went wrong different size of time series.";
//    std::exit(EXIT_FAILURE);
//  }
}
//
//data_HB_1d::data_HB_1d(unsigned _ndata)
//{
//
//  Prec.resize(_ndata, 999999);
//  AEt = Prec;
//  PEt = Prec;
//  Temp = Prec;
//  SoiS = Prec;
//  GroS = Prec;
//  SurS = Prec;
//  IntS = Prec;
//  TotR = Prec;
//  Basf = Prec;
//  DirR = Prec;
//  OveF = Prec;
//  Infi = Prec;
//  Perc = Prec;
//  Pref = Prec;
//
//  std::cout << "\nSize of Perc : " << Perc.size() << "\n";
//  std::cout << "\nSize of Perc : " << Perc.size() << " " << Prec[_ndata-1] << "\n";
//
//
//  if(Prec.size()!= Perc.size()) {
//        std::cout << "\n Something went wrong different size of time series.";
//        std::exit(EXIT_FAILURE);
//  }
//
//  numTS = Pref.size();
//
//}
//
//
data_HB_1d::~data_HB_1d() {
  //dtor

}

data_HB_1d::data_HB_1d(const data_HB_1d& other): numTS(0),
  Latitude(0.0),
  PETtype(pet_Type::OUDIN),
  year(1,1),
  month(1,1),
  day(1,1),
  Jday(1,1),
  init_year(1,1),
  init_month(1,1),
  init_day(1,1),
  Prec(1,1),
  Snow(1,1),
  AEt(1,1),
  PEt(1,1),
  Temp(1,1),
  TroF(1,1),
  SteF(1,1),
  CanF(1,1),
  CanS(1,1),
  SteS(1,1),
  EvaC(1,1),
  EvaS(1,1),
  EvbS(1,1),
  IntS(1,1),
  SoiS(1,1),
  GroS(1,1),
  SurS(1,1),
  TotR(1,1),
  Basf(1,1),
  DirR(1,1),
  Melt(1,1),
  Perc(1,1),
  Pref(1,1),
  Etsw(1,1),
  PonS(1,1),
  EtpO(1,1),
  PoiS(1,1),
  PoiG(1,1),
  init_SoiS(0.0),
  init_GroS(0.0),
  init_CanS(0.0),
  init_SteS(0.0),
  init_SnoS(0.0),
  init_SurS(0.0),
  init_GroS1(0.0),
  init_GroS2(0.0),
  numfastRes(1),
  StateFastRes(1,1),
  OutFastRes(1,1) {
  numTS = other.numTS;
  Latitude = other.Latitude;
  PETtype = other.PETtype;
  year = other.year;
  month = other.month;
  day = other.day;
  Jday = other.Jday;
  init_year = other.init_year;
  init_month = other.init_month;
  init_day = other.init_day;
  Prec = other.Prec;//!< Precipitation
  Snow = other.Snow;//!< Snow depth
  AEt = other.AEt;//!< Actual Evapotranspiration
  PEt = other.PEt;//!< Potential Evapotranspiration
  Temp = other.Temp;//!< Temperature
  TroF = other.TroF;//!< Throughfall
  SteF = other.SteF;//!< Stemflow
  CanF = other.CanF;//!< Canopy drainage
  CanS = other.CanS;//!< Canopy storage
  SteS = other.SteS;//!< Stem storage
  EvaC = other.EvaC;//!< Canopy Evaporation
  EvaS = other.EvaS;//!< Stem Evaporation
  EvbS = other.EvbS;//!< Bare soil Evapotranspiration
  IntS = other.IntS;//!< Interception storage
  SoiS = other.SoiS;//!< Soil storage
  GroS = other.GroS;//!< Groundwater storage
  SurS = other.SurS;//!< Surface retention
  TotR = other.TotR;//!< Total runoff
  Basf = other.Basf;//!< Baseflow
  DirR = other.DirR;//!< Direct Runoff
  Melt = other.Melt;//!< Melting
  Perc = other.Perc;//!< Percolation
  Pref = other.Pref;//!< Effective Precipitation
  Etsw = other.Etsw;//!< Evaporation from surface retention
  EtpO = other.EtpO;//!< Evaporation from pond surface [mm/day]
  PonS = other.PonS;//!< Pond storage
  PoiS = other.PoiS;//!< Percolation between the pond and the soil
  PoiG = other.PoiG;//!< Percolation between the pond and the groundwater
  init_SoiS = other.init_SoiS;//!< Initial value of soil storage
  init_GroS = other.init_GroS;//!< Initial value of groundwater storage
  init_CanS = other.init_SteS;//!< Initial value of Canopy Interception storage
  init_SteS = other.init_SteS;//!< Initial value of Stem Interception storage
  init_SnoS = other.init_SnoS;//!< Initial variable of Snow storage
  init_SurS = other.init_SurS;//!< Initial value of Surface retention storage
  init_GroS1 = other.init_GroS1;
  init_GroS2 = other.init_GroS2;
  numfastRes = other.numfastRes;
  StateFastRes = other.StateFastRes;
  OutFastRes = other.OutFastRes;
}

data_HB_1d& data_HB_1d::operator=(const data_HB_1d& rhs) {

  if (this == &rhs) return *this;
  else {
    numTS = rhs.numTS;
    Latitude = rhs.Latitude;
    PETtype = rhs.PETtype;
    year = rhs.year;
    month = rhs.month;
    day = rhs.day;
    Jday = rhs.Jday;
    init_year = rhs.init_year;
    init_month = rhs.init_month;
    init_day = rhs.init_day;
    Prec = rhs.Prec;//!< Precipitation
    Snow = rhs.Snow;//!< Snow depth
    AEt = rhs.AEt;//!< Actual Evapotranspiration
    PEt = rhs.PEt;//!< Potential Evapotranspiration
    Temp = rhs.Temp;//!< Temperature
    TroF = rhs.TroF;//!< Throughfall
    SteF = rhs.SteF;//!< Stemflow
    CanF = rhs.CanF;//!< Canopy drainage
    CanS = rhs.CanS;//!< Canopy storage
    SteS = rhs.SteS;//!< Stem storage
    EvaC = rhs.EvaC;//!< Canopy Evaporation
    EvaS = rhs.EvaS;//!< Stem Evaporation
    EvbS = rhs.EvbS;//!< Bare soil Evapotranspiration
    IntS = rhs.IntS;//!< Interception storage
    SoiS = rhs.SoiS;//!< Soil storage
    GroS = rhs.GroS;//!< Groundwater storage
    SurS = rhs.SurS;//!< Surface retention
    TotR = rhs.TotR;//!< Total runoff
    Basf = rhs.Basf;//!< Baseflow
    DirR = rhs.DirR;//!< Direct Runoff
    Melt = rhs.Melt;//!< Melting
    Perc = rhs.Perc;//!< Percolation
    Pref = rhs.Pref;//!< Effective Precipitation
    Etsw = rhs.Etsw;//!< Evaporation from surface retention
    EtpO = rhs.EtpO;//!< Evaporation from pond surface [mm/day]
    PonS = rhs.PonS;//!< Pond storage
    PoiS = rhs.PoiS;//!< Percolation between the pond and the soil
    PoiG = rhs.PoiG;//!< Percolation between the pond and the groundwater
    init_SoiS = rhs.init_SoiS;//!< Initial value of soil storage
    init_GroS = rhs.init_GroS;//!< Initial value of groundwater storage
    init_CanS = rhs.init_SteS;//!< Initial value of Canopy Interception storage
    init_SteS = rhs.init_SteS;//!< Initial value of Stem Interception storage
    init_SnoS = rhs.init_SnoS;//!< Initial variable of Snow storage
    init_SurS = rhs.init_SurS;//!< Initial value of Surface retention storage
    init_GroS1 = rhs.init_GroS1;
    init_GroS2 = rhs.init_GroS2;
    numfastRes = rhs.numfastRes;
    StateFastRes = rhs.StateFastRes;
    OutFastRes = rhs.OutFastRes;
  }
  return *this;
}

void data_HB_1d::s_data(const hdata& dta,const ts_type& _tsType, bool updateNumTS) {

  if(updateNumTS)  numTS = dta.size();

  switch (_tsType) {
  case ts_type::PREC:
    Prec = dta;
//    std::cout << "New precipitation --> loaded\n";
    break;
  case ts_type::SNOW:
    Snow = dta;
//    std::cout << "New snow --> loaded\n";
    break;
  case ts_type::AET:
    AEt = dta;
//    std::cout << "New AET --> loaded\n";
    break;
  case ts_type::PET:
    PEt = dta;
//    std::cout << "New PET --> loaded\n";
    break;
  case ts_type::TEMP:
    Temp = dta;
//    std::cout << "New TEMP --> loaded\n";
    break;
  case ts_type::TROF:
    TroF = dta;
//    std::cout << "New TROF --> loaded\n";
    break;
  case ts_type::STEF:
    SteF = dta;
//    std::cout << "New STEF --> loaded\n";
    break;
  case ts_type::CANF:
    CanF = dta;
//    std::cout << "New CANF --> loaded\n";
    break;
  case ts_type::CANS:
    CanS = dta;
//    std::cout << "New CANS --> loaded\n";
    break;
  case ts_type::STES:
    SteS = dta;
//    std::cout << "New STES --> loaded\n";
    break;
  case ts_type::EVAC:
    EvaC = dta;
//    std::cout << "New EVAC --> loaded\n";
    break;
  case ts_type::EVAS:
    EvaS = dta;
//    std::cout << "New EVAS --> loaded\n";
    break;
  case ts_type::EVBS:
    EvbS = dta;
//    std::cout << "New EVBS --> loaded\n";
    break;
  case ts_type::SOIS:
    SoiS = dta;
//    std::cout << "New SoiS --> loaded\n";
    break;
  case ts_type::GROS:
    GroS = dta;
//    std::cout << "New GroS --> loaded\n";
    break;
  case ts_type::SURS:
    SurS= dta;
//    std::cout << "New SurS--> loaded\n";
    break;
  case ts_type::INTS:
    IntS = dta;
//    std::cout << "New IntS --> loaded\n";
    break;
  case ts_type::TOTR:
    TotR = dta;
//    std::cout << "New TotR --> loaded\n";
    break;
  case ts_type::BASF:
    Basf = dta;
//    std::cout << "New Basf --> loaded\n";
    break;
  case ts_type::DIRR:
    DirR = dta;
//    std::cout << "New DirR --> loaded\n";
    break;
  case ts_type::MELT:
    Melt = dta;
//    std::cout << "New Melt --> loaded\n";
    break;
  case ts_type::PERC:
    Perc = dta;
//    std::cout << "New Perc --> loaded\n";
    break;
  case ts_type::PREF:
    Pref = dta;
//    std::cout << "New Peef --> loaded\n";
    break;
  case ts_type::ETSW:
    Etsw = dta;
    //    std::cout << "New ETSW --> loaded\n";
    break;
  case ts_type::PONS:
    PonS = dta;
    //    std::cout << "New PONS --> loaded\n";
    break;
  case ts_type::ETPO:
    EtpO = dta;
    //    std::cout << "New ETPO --> loaded\n";
    break;
  case ts_type::POIS:
    PoiS = dta;
    //    std::cout << "New POIS--> loaded\n";
    break;
  case ts_type::POIG:
    PoiG = dta;
    //    std::cout << "New POIG --> loaded\n";
    break;
  }
//    //PRE.insert(PRE.begin(), data);
//    if(!Erase) {
//        PRE.at(pos) =  data;
//      } else {
//          PRE.clear();
//
//      }
//
//
//    std::cout  << std::endl << PRE.size() << "\t" << PRE[0] << std::endl;
}

void data_HB_1d::s_varVal(const numberSel& dta, const unsigned& tst,const ts_type& _tsType) {

  switch (_tsType) {
  case ts_type::PREC:
    Prec[tst] = dta;
    break;
  case ts_type::SNOW:
    Snow[tst] = dta;
    break;
  case ts_type::AET:
    AEt[tst] = dta;
    break;
  case ts_type::PET:
    PEt[tst] = dta;
    break;
  case ts_type::TEMP:
    Temp[tst] = dta;
    break;
  case ts_type::TROF:
    TroF[tst] = dta;
    break;
  case ts_type::STEF:
    SteF[tst] = dta;
    break;
  case ts_type::CANF:
    CanF[tst] = dta;
    break;
  case ts_type::CANS:
    CanS[tst] = dta;
    break;
  case ts_type::STES:
    SteS[tst] = dta;
    break;
  case ts_type::EVAC:
    EvaC[tst] = dta;
    break;
  case ts_type::EVAS:
    EvaS[tst] = dta;
    break;
  case ts_type::EVBS:
    EvbS[tst] = dta;
    break;
  case ts_type::SOIS:
    SoiS[tst] = dta;
    break;
  case ts_type::GROS:
    GroS[tst] = dta;
    break;
  case ts_type::SURS:
    SurS[tst] = dta;
    break;
  case ts_type::INTS:
    IntS[tst] = dta;
    break;
  case ts_type::TOTR:
    TotR[tst] = dta;
    break;
  case ts_type::BASF:
    Basf[tst] = dta;
    break;
  case ts_type::DIRR:
    DirR[tst] = dta;
    break;
  case ts_type::MELT:
    Melt[tst] = dta;
    break;
  case ts_type::PERC:
    Perc[tst] = dta;
    break;
  case ts_type::PREF:
    Pref[tst] = dta;
    break;
  case ts_type::ETSW:
    Etsw[tst] = dta;
    break;
  case ts_type::PONS:
    PonS[tst] = dta;
    break;
  case ts_type::ETPO:
    EtpO[tst] = dta;
    break;
  case ts_type::POIS:
    PoiS[tst] = dta;
    break;
  case ts_type::POIG:
    PoiG[tst] = dta;
    break;

}
  return ;


}


numberSel data_HB_1d::g_dta(const unsigned& tst,const ts_type& _tsType) {

  switch (_tsType) {
  case ts_type::PREC:
    return Prec[tst];
  case ts_type::SNOW:
    return Snow[tst];
  case ts_type::AET:
    return AEt[tst];
  case ts_type::PET:
    return PEt[tst];
  case ts_type::TEMP:
    return Temp[tst];
  case ts_type::TROF:
    return TroF[tst];
  case ts_type::STEF:
    return SteF[tst];
  case ts_type::CANF:
    return CanF[tst];
  case ts_type::EVAC:
    return EvaC[tst];
  case ts_type::EVAS:
    return EvaS[tst];
  case ts_type::EVBS:
    return EvbS[tst];
  case ts_type::CANS:
    return CanS[tst];
  case ts_type::STES:
    return SteS[tst];
  case ts_type::SOIS:
    return SoiS[tst];
  case ts_type::GROS:
    return GroS[tst];
  case ts_type::SURS:
    return SurS[tst];
  case ts_type::INTS:
    return IntS[tst];
  case ts_type::TOTR:
    return TotR[tst];
  case ts_type::BASF:
    return Basf[tst];
  case ts_type::DIRR:
    return DirR[tst];
  case ts_type::MELT:
    return Melt[tst];
  case ts_type::PERC:
    return Perc[tst];
  case ts_type::PREF:
    return Pref[tst];
  case ts_type::ETSW:
    return Etsw[tst];
  case ts_type::PONS:
    return PonS[tst];
   case ts_type::ETPO:
    return EtpO[tst];
  case ts_type::POIS:
    return PoiS[tst];
  case ts_type::POIG:
    return PoiG[tst];
  }

  return 0;
}


void data_HB_1d::s_initStates(const hdata& initfastRes, const numberSel& init_State,const init_Stype& _Stype) {

  switch (_Stype) {
  case init_Stype::SOIL:
    init_SoiS = init_State;
    break;
  case init_Stype::GROUNDWAT:
    init_GroS = init_State;
    break;
  case init_Stype::CANS:
    init_CanS = init_State;
    break;
  case init_Stype::STES:
    init_SteS = init_State;
    break;
  case init_Stype::SNOS:
    init_SnoS = init_State;
    break;
  case init_Stype::SURFRET:
    init_SurS = init_State;
    break;
  case init_Stype::FASTRES:
    s_init_states_fastRes(initfastRes);
    break;
  case init_Stype::GROS1:
    init_GroS1 = init_State;
    break;
  case init_Stype::GROS2:
    init_GroS2 = init_State;
    break;
  }

  return ;
}


void data_HB_1d::s_init_states_fastRes(const hdata& initfastRes) {

  StateFastRes = initfastRes;
  OutFastRes = initfastRes;

}

numberSel data_HB_1d::g_initState(const init_Stype& _Stype) {

  switch (_Stype) {
  case init_Stype::SOIL:
    return init_SoiS;
    break;
  case init_Stype::GROUNDWAT:
    return init_GroS;
    break;
  case init_Stype::CANS:
    return init_CanS;
    break;
  case init_Stype::STES:
    return init_SteS;
    break;
  case init_Stype::SNOS:
    return init_SnoS;
    break;
  case init_Stype::SURFRET:
    return init_SurS;
    break;
  case init_Stype::GROS1:
    return init_GroS1;
    break;
  case init_Stype::GROS2:
    return init_GroS2;
    break;
  case init_Stype::FASTRES:
    return StateFastRes[1];// returns only the init state of the first reservoir
    break;
  }

  return 0;
}

hdata data_HB_1d::g_initFastresStates() {

  return StateFastRes;

}


numberSel data_HB_1d::g_oneFastResOut(const unsigned& itFasRes) {

  return OutFastRes[itFasRes];

}

numberSel data_HB_1d::g_oneFastResState(const unsigned& itFasRes) {

  return StateFastRes[itFasRes];

}


void data_HB_1d::s_oneFastreResUout(const numberSel& resOut,const unsigned& itFasRes) {

  OutFastRes[itFasRes] = resOut;

  return ;
}


void data_HB_1d::s_oneFastreResUstate(const numberSel& resOut,const unsigned& itFasRes) {

  StateFastRes[itFasRes] = resOut;

  return ;
}


unsigned data_HB_1d::g_numdta() {

  return numTS;

}


void data_HB_1d::calc_Pet() {

  switch (PETtype) {
  case pet_Type::OUDIN:
    OudinPET();
    break;
  case pet_Type::HAMON:
    HamonPET();
    break;
  case pet_Type::THORNTHWAITE:
    ThornthwaitePET();
    break;
  case pet_Type::BLANEYCRIDDLE:
    BlaneycriddlePET();
    break;
  case pet_Type::JENSENHAISE:
    JensenhaisePET();
    break;
  case pet_Type::MCGUINNESSBORDNE:
    McGuinnessbrdnePET();
    break;
  }
  return ;
}



void data_HB_1d::s_Pet_Pars(const numberSel& newLatitude, const pet_Type& newPeType) {

  Latitude = newLatitude;
  PETtype = newPeType;

//  s_Julian_day();

//  calc_Pet();

  return ;
}


void data_HB_1d::OudinPET() {

  numberSel dec = 0.0, rad_Latid = Latitude / 180 * M_PI, ndy = 365, dr = 0.0, omega = 0.0, Ra = 0.0, GlSoCo = 0.0820;

  for(unsigned tst=0; tst<numTS; tst++) {
    if( leap_Check_Year(year[tst]) ) ndy = 366;
    else ndy = 365;
    dec = 0.409 * sin(2 * (static_cast<numberSel>(Jday[tst])) * M_PI / ndy - 1.39);
    dr = 1 + 0.033 * cos((static_cast<numberSel>(Jday[tst])) * 2 * M_PI / ndy);
    omega = acos(-tan(rad_Latid) * tan(dec));
    Ra = (24 * 60) / M_PI * GlSoCo * dr * (omega * sin(rad_Latid) * sin(dec) + cos(rad_Latid) * cos(dec) * sin(omega));
    if ((Temp[tst] + 5.0) >= 0.0 ) PEt[tst] = 0.408 * Ra * (Temp[tst] + 5) / 100;
    else PEt[tst] = 0.0;
  }

  PETtype = pet_Type::OUDIN;

  return ;
}

/** \brief The implementation of Hamon calculation for PET according to
 *
 *
 */

void data_HB_1d::HamonPET() {

  numberSel dec = 0.0, rad_Latid = Latitude / 180 * M_PI, ndy = 365, Nn = 0.0, omega = 0.0, Esat = 0.0;

  for(unsigned tst=0; tst<numTS; tst++) {
    if(leap_Check_Year(year[tst])) ndy = 366;
    else ndy = 365;
    dec = 0.409 * sin( (2 * M_PI ) * static_cast<numberSel>(Jday[tst]) / ndy - 1.39);
    omega = acos( (-1)* tan(rad_Latid) * tan(dec));
    Nn = 24 / M_PI * omega;
    Esat = 6.108 * exp( 17.27 * Temp[tst] / (Temp[tst]+237.3) );
    PEt[tst] = 0.1651 * 216.7 * (Nn/12) * ( Esat / (Temp[tst]+273.3) );
  }

  PETtype = pet_Type::HAMON;

  return ;

}


void data_HB_1d::ThornthwaitePET() {

  numberSel rad_Latid = Latitude / 180 * M_PI, ndy = 365, omega = 0.0, dec = 0.0;
  hdata Nn(1.01,1);
  Nn.resize(numTS);

  for(unsigned tst=0; tst<numTS; tst++) {
    if(leap_Check_Year(year[tst])) ndy = 366;
    else ndy = 365;
    dec = 0.409 * sin( (2 * M_PI ) * static_cast<numberSel>(Jday[tst]) / ndy - 1.39);
    omega = acos( (-1)* tan(rad_Latid) * tan(dec));
    Nn[tst] = 24 / M_PI * omega;
  }

  unsigned nmonthsinData = 1, nyearsinData = 1;

  for(unsigned tst=1; tst<numTS; tst++) {
    if(month[tst] != month[tst-1]) nmonthsinData++;
    if(year[tst] != year[tst-1]) nyearsinData++;
  }

  nmonthsinData++;//last month must be considered
  nyearsinData++;//last year must be considered

  hdata tam(1.01,1), Nmeanmonth(1.01,1), numDaysMonth(1.01,1);//monthly temperatures
  tam.resize(nmonthsinData);
  Nmeanmonth.resize(nmonthsinData);
  numDaysMonth.resize(nmonthsinData);
  hdata helpyear(1.01,1);
  helpyear.resize(nmonthsinData);

  numberSel helptam =0, helpNn=0, helpNdaysInMoth = 1;
  helptam = helptam + Temp[0];
  helpNn = helpNn + Nn[0];

  unsigned helpInd =0;
  for(unsigned tst=1; tst<numTS; tst++) {
    if((month[tst] > month[tst-1])||(year[tst] > year[tst-1])) {
      helptam = 0.0;
      helpNn =  0.0;
      helpInd++;
      helpNdaysInMoth = 0;
      }
    helptam = helptam + Temp[tst];
    helpNn = helpNn + Nn[tst];
    helpNdaysInMoth = helpNdaysInMoth + 1.0;
    tam[helpInd] = helptam / helpNdaysInMoth;
    Nmeanmonth[helpInd] = helpNn / helpNdaysInMoth;
    numDaysMonth[helpInd] = helpNdaysInMoth;
    helpyear[helpInd] = (numberSel) year[tst];
  }

  hdata i_heatindex(1.01,1);
  i_heatindex.resize(nmonthsinData);
  for(unsigned it=0; it<nmonthsinData;it++){
    if(tam[it]<0.0) {
      tam[it] = 0.0;
    } else {
      if(tam[it]>0) {
        i_heatindex[it] = std::pow((tam[it]*0.2),1.514);
      } else i_heatindex[it]=0.0;
      }
  }

  hdata I_annualHeatindex(1.01,1);
  I_annualHeatindex.resize(nyearsinData);
  hdata acoeff(1.01,1);
  acoeff.resize(nyearsinData);

  numberSel helpsumI = 0.0, part1 = 0.0, part2 =0.0, part0 = 0.0;
  unsigned helpIn = 0;

  numberSel lastyear = helpyear.max();
  helpsumI = i_heatindex[0];
  for(unsigned it=1; it<nmonthsinData; it++){
    helpsumI = helpsumI + i_heatindex[it];
    I_annualHeatindex[helpIn] = helpsumI - i_heatindex[it];
    part1 = I_annualHeatindex[helpIn] * I_annualHeatindex[helpIn]* I_annualHeatindex[helpIn];
    part2 = I_annualHeatindex[helpIn] * I_annualHeatindex[helpIn];
    part0 = 0.01729 * I_annualHeatindex[helpIn];
    acoeff[helpIn] = 0.0000006751 * part1 - 0.0000771 * part2 +  part0 + 0.49239;
    if((helpyear[it] != helpyear[it-1])){
      helpIn++;
      helpsumI = i_heatindex[it];
      }
    }

  hdata amonthly(1.01,1);
  hdata Imonthly(1.01,1);
  hdata Epetraw(0.01,1);
  amonthly.resize(nmonthsinData);
  Imonthly.resize(nmonthsinData);
  Epetraw.resize(nmonthsinData);
  amonthly[0] = acoeff[0];
  Imonthly[0] = I_annualHeatindex[0];
  Epetraw[0] = 16*(pow((10*tam[0]/Imonthly[0]),amonthly[0]));

  unsigned helpYearInd = 0.0;
  for(unsigned it=1; it<nmonthsinData; it++){
    if(helpyear[it] != helpyear[it-1]){
      helpYearInd++;
      }
    amonthly[it] = acoeff[helpYearInd];
    Imonthly[it] = I_annualHeatindex[helpYearInd];
    Epetraw[it] = 16*(pow((10*tam[it]/Imonthly[it]),amonthly[it]));
  }

  hdata Epet(1.01,1);
  Epet.resize(nmonthsinData);
  for(unsigned it=0; it<nmonthsinData; it++){
    Epet[it] = Epetraw[it] * (Nmeanmonth[it] / 12) * (numDaysMonth[it] / 30);
  }

  unsigned helpit =0;
  PEt[0] = Epet[0] / numDaysMonth[0];
  // numberSel difPETm=0.0, cdifPet = 0.0;
  // difPETm = Epet[0] - Epet[1] / numDaysMonth[0];
  // cdifPet = 0.0;
  for(unsigned tst=1; tst<numTS; tst++) {
    if(month[tst]!=month[tst-1]){
      helpit++;
      // difPETm =(Epet[helpit] - Epet[helpit+1]) / numDaysMonth[helpit];
      // cdifPet = 0.0;
      }
    // cdifPet = cdifPet + difPETm;
    PEt[tst] = Epet[helpit] / numDaysMonth[helpit];
    // PEt[tst] = Epet[helpit] + cdifPet;
  }

  PETtype = pet_Type::THORNTHWAITE;

  return ;
}

void data_HB_1d::BlaneycriddlePET(){

  numberSel rad_Latid = Latitude / 180 * M_PI, ndy = 365, omega = 0.0, dec = 0.0;
  hdata Nn(1.01,1);
  Nn.resize(numTS);

  for(unsigned tst=0; tst<numTS; tst++) {
    if(leap_Check_Year(year[tst])) ndy = 366;
    else ndy = 365;
    dec = 0.409 * sin( (2 * M_PI ) * static_cast<numberSel>(Jday[tst]) / ndy - 1.39);
    omega = acos( (-1)* tan(rad_Latid) * tan(dec));
    Nn[tst] = 24 / M_PI * omega / (365*12);
  }

  for(unsigned tst=1; tst<numTS; tst++) {
    PEt[tst] = (Nn[tst] * 0.85 ) * 100 * (0.46 * Temp[tst] + 8.13);
    //Xu and Singh 2001
    //Evaluation and generalization of temperature-basedmethods for
    //calculating evaporation Hydrol. Process. 15, 305–319 (2001)
    }

  return ;
}


void data_HB_1d::JensenhaisePET(){
  numberSel dec = 0.0, rad_Latid = Latitude / 180 * M_PI, ndy = 365, dr = 0.0, omega = 0.0, Ra = 0.0, GlSoCo = 0.0820;

  for(unsigned tst=0; tst<numTS; tst++) {
    if( leap_Check_Year(year[tst]) ) ndy = 366;
    else ndy = 365;
    dec = 0.409 * sin(2 * (static_cast<numberSel>(Jday[tst])) * M_PI / ndy - 1.39);
    dr = 1 + 0.033 * cos((static_cast<numberSel>(Jday[tst])) * 2 * M_PI / ndy);
    omega = acos(-tan(rad_Latid) * tan(dec));
    Ra = (24 * 60) / M_PI * GlSoCo * dr * (omega * sin(rad_Latid) * sin(dec) + cos(rad_Latid) * cos(dec) * sin(omega));
    if (1000 * Ra * Temp[tst] / (40 * 2450) <0) PEt[tst] =0.0;
      else PEt[tst] = 1000 * Ra * Temp[tst] / (40 * 2450);
  }
  //Oudin 2005 Joh

  return ;
};


void data_HB_1d::McGuinnessbrdnePET(){
  numberSel dec = 0.0, rad_Latid = Latitude / 180 * M_PI, ndy = 365, dr = 0.0, omega = 0.0, Ra = 0.0, GlSoCo = 0.0820;

  for(unsigned tst=0; tst<numTS; tst++) {
    if( leap_Check_Year(year[tst]) ) ndy = 366;
    else ndy = 365;
    dec = 0.409 * sin(2 * (static_cast<numberSel>(Jday[tst])) * M_PI / ndy - 1.39);
    dr = 1 + 0.033 * cos((static_cast<numberSel>(Jday[tst])) * 2 * M_PI / ndy);
    omega = acos(-tan(rad_Latid) * tan(dec));
    Ra = (24 * 60) / M_PI * GlSoCo * dr * (omega * sin(rad_Latid) * sin(dec) + cos(rad_Latid) * cos(dec) * sin(omega));
    if (1000 * Ra * (Temp[tst] +5) / (68 * 2450) <0) PEt[tst] =0.0;
    else PEt[tst] = 1000 * Ra * (Temp[tst] +5) / (68 * 2450);
  }
  //Oudin 2005 Joh

  return ;
}
/** \brief The description  of leap years
 *
 * \param year to be described
 *
 */
bool data_HB_1d::leap_Check_Year(unsigned TestedYear) {

  return (((TestedYear % 4 == 0) && (TestedYear % 100 != 0)) || (TestedYear % 400 == 0));

}

/** \brief Setting the rank of a day (almost as Julian day not translated by CE constant)
 *
 *  Ranks the day for pet calculations
 */
void data_HB_1d::s_Julian_day() {

  caldata cumNdayinYear{0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  caldata cumNdayinYearLeafyear {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};

  for( unsigned tst=0; tst<numTS; tst++) {
    if (leap_Check_Year(year[tst])) {
      Jday[tst] = cumNdayinYearLeafyear[ month[tst]-1] + day[tst];
    } else {
      Jday[tst] = cumNdayinYear[ month[tst]-1] + day[tst];
    }
  }
}

/** \brief Setting the whole calender
 *
 *
 */
void data_HB_1d::s_calender() {

  caldata NdaysInYear{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  caldata NdaysInYearLeapYear {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  year.resize(numTS);
  month.resize(numTS);
  day.resize(numTS);
  Jday.resize(numTS);

  year[0] = init_year[0];
  month[0] = init_month[0];
  day[0] = init_day[0];
  for(unsigned dt=1; dt<numTS; dt++) {
    if((month[dt-1] ==12) && (day[dt-1] == 31)) {
      day[dt] = 1;
      month[dt] = 1;
      year[dt] = year[dt-1] + 1;
    } else {
      if(leap_Check_Year(year[dt-1])) {
        if(day[dt-1] == NdaysInYearLeapYear[month[dt-1]-1]) {
          day[dt] = 1;
          month[dt] = month[dt-1] + 1;
          year[dt] = year[dt-1];
        } else {
          day[dt] = day[dt-1] +1;
          month[dt] = month[dt-1];
          year[dt] = year[dt-1];
        }
      } else {
        if(day[dt-1] == NdaysInYear[month[dt-1]-1]) {
          day[dt] = 1;
          month[dt] = month[dt-1] + 1;
          year[dt] = year[dt-1];
        } else {
          day[dt] = day[dt-1] +1;
          month[dt] = month[dt-1];
          year[dt] = year[dt-1];
        }
      }
    }
  }

  s_Julian_day();

  return ;
}

numberSel data_HB_1d::get_daysInMonth(const unsigned& tstMonth, const unsigned& year){

  numberSel ndays = 0;
  hdata NdaysInYear{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  hdata NdaysInYearLeapYear {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};


  if(leap_Check_Year(year)) {
    ndays = NdaysInYearLeapYear[tstMonth-1];
    }
  else {
    ndays = NdaysInYear[tstMonth-1];
    }


  return ndays;
}


/** \brief Printing the calender of data
 *
 *
 */
// void data_HB_1d::p_calender() {
//
//   std::cout << std::endl << "Printing the calender" << std::endl << std::endl;
//   for (unsigned dt=0; dt<numTS ; dt++ )
//     std::cout << year[dt] << "\t" << month[dt] << "\t" << day[dt] << "\t" << Jday[dt] << std:: endl;
//
// }

/** \brief Initialization of the first Date
 *
 *
 */
void data_HB_1d::s_initDate(const unsigned& Year, const unsigned& Month,const unsigned& Day,const unsigned& initNumTS) {

  init_year = Year;
  init_month = Month;
  init_day = Day;
  numTS = initNumTS;

}


numberDta data_HB_1d::g_calDta(const cal_Type& calDate, const unsigned& ts) {

  switch (calDate) {
  case cal_Type::YEAR:
    return year[ts];
    break;
  case cal_Type::MONTH:
    return month[ts];
    break;
  case cal_Type::DAY:
    return day[ts];
    break;
  case cal_Type::JDAY:
    return Jday[ts];
    break;
  }

  return 0;

}

hdata data_HB_1d::get_HbTsData(const ts_type& _tsType) {

  switch (_tsType) {
  case ts_type::PREC:
    return Prec;
  case ts_type::SNOW:
    return Snow;
  case ts_type::AET:
    return AEt;
  case ts_type::PET:
    return PEt;
  case ts_type::TEMP:
    return Temp;
  case ts_type::TROF:
    return TroF;
  case ts_type::STEF:
    return SteF;
  case ts_type::CANF:
    return CanF;
  case ts_type::EVAC:
    return EvaC;
  case ts_type::EVAS:
    return EvaS;
  case ts_type::EVBS:
    return EvbS;
  case ts_type::CANS:
    return CanS;
  case ts_type::STES:
    return SteS;
  case ts_type::SOIS:
    return SoiS;
  case ts_type::GROS:
    return GroS;
  case ts_type::SURS:
    return SurS;
  case ts_type::INTS:
    return IntS;
  case ts_type::TOTR:
    return TotR;
  case ts_type::BASF:
    return Basf;
  case ts_type::DIRR:
    return DirR;
  case ts_type::MELT:
    return Melt;
  case ts_type::PERC:
    return Perc;
  case ts_type::PREF:
    return Pref;
  case ts_type::ETSW:
    return Etsw;
  case ts_type::PONS:
    return PonS;
  case ts_type::ETPO:
    return EtpO;
  case ts_type::POIS:
    return PoiS;
  case ts_type::POIG:
    return PoiG;
  }

  hdata helpV(-9999,1);

  return helpV;

}


void data_HB_1d::setAllToZeros(const bool& CalDta,const bool& TsDta) {

  if(CalDta) {
    for(const auto& it : all_caDT) {
      setOneCalDateToZero(it);
    }
  }

  if(TsDta) {
    for(const auto& it : all_ts) {
      setOneTstoZero(it);
    }
  }

  return ;

}


void data_HB_1d::setOneTstoZero(const ts_type& _tsType) {

  switch (_tsType) {
  case ts_type::PREC:
    Prec = 0.0;
    break;
  case ts_type::SNOW:
    Snow = 0.0;
    break;
  case ts_type::AET:
    AEt = 0.0;
    break;
  case ts_type::PET:
    PEt = 0.0;
    break;
  case ts_type::TEMP:
    Temp = 0.0;
    break;
  case ts_type::TROF:
    TroF = 0.0;
    break;
  case ts_type::STEF:
    SteF = 0.0;
    break;
  case ts_type::CANF:
    CanF = 0.0;
    break;
  case ts_type::CANS:
    CanS = 0.0;
    break;
  case ts_type::STES:
    SteS = 0.0;
    break;
  case ts_type::EVAC:
    EvaC = 0.0;
    break;
  case ts_type::EVAS:
    EvaS = 0.0;
    break;
  case ts_type::EVBS:
    EvbS = 0.0;
    break;
  case ts_type::SOIS:
    SoiS = 0.0;
    break;
  case ts_type::GROS:
    GroS = 0.0;
    break;
  case ts_type::SURS:
    SurS = 0.0;
    break;
  case ts_type::INTS:
    IntS = 0.0;
    break;
  case ts_type::TOTR:
    TotR = 0.0;
    break;
  case ts_type::BASF:
    Basf = 0.0;
    break;
  case ts_type::DIRR:
    DirR = 0.0;
    break;
  case ts_type::MELT:
    Melt = 0.0;
    break;
  case ts_type::PERC:
    Perc = 0.0;
    break;
  case ts_type::PREF:
    Pref = 0.0;
    break;
  case ts_type::ETSW:
    Etsw = 0.0;
    break;
  case ts_type::PONS:
    PonS = 0.0;
    break;
  case ts_type::ETPO:
    EtpO = 0.0;
    break;
  case ts_type::POIS:
    PoiS = 0.0;
    break;
  case ts_type::POIG:
    PoiG = 0.0;
    break;
  }

  return ;

}


void data_HB_1d::setOneCalDateToZero(const cal_Type& calDate) {

  switch (calDate) {
  case cal_Type::YEAR:
    year = 0;
    break;
  case cal_Type::MONTH:
    month = 0;
    break;
  case cal_Type::DAY:
    day = 0;
    break;
  case cal_Type::JDAY:
    Jday = 0;
    break;
  }

  return ;

}

void data_HB_1d::printDataToFile(const std::string& Filet) {

  std::ofstream outfilet;
  outfilet.open (Filet.c_str());

  if(outfilet.is_open()) {
    outfilet.precision(6);

    outfilet << std::right;
    char setfiller = ' ';
    for(unsigned Ts=0; Ts<numTS; Ts++) {
      for(const auto& it : all_caDT) {
        outfilet << std::setprecision(0);
        outfilet << std::setw(2);
        outfilet << g_calDta(it,Ts) << " ";
//        std::cout << g_calDta(it,Ts) << " ";
      }
//     std::cout  << std::endl;
      for(const auto& it : all_ts) {
        outfilet << std::setprecision(2);
        outfilet << std::fixed << std::setw(10);
        outfilet << std::setfill(setfiller) << std::right;
        outfilet << g_dta(Ts, it) << " ";
      }
      outfilet << std::endl;
    }
    outfilet.close();

  }
  // else {
  //
  //   std::cout << "There is a error in opening and writing to the file with path and name: " << Filet.c_str() << std::endl;
  //
  // }


}

caldata data_HB_1d::getCalData(const cal_Type& calDate){

  switch (calDate) {
  case cal_Type::YEAR:
    return year;
    break;
  case cal_Type::MONTH:
    return month;
    break;
  case cal_Type::DAY:
    return day;
    break;
  case cal_Type::JDAY:
    return Jday;
    break;
  }

  caldata helpV(-9999,1);

  return helpV;
}


void data_HB_1d::loadCalData(const caldata& yyear, const caldata& mmonth, const caldata& dday) {

  year = yyear;
  month = mmonth;
  day = dday;

  Jday = day;

  s_Julian_day();

  return ;

}

void data_HB_1d::s_numFastRes(const numberDta& nFastRes) {

  numfastRes = nFastRes;

  return ;
}

numberDta data_HB_1d::g_numFastRes() {

  return numfastRes;

}
