#include "dta_dam_1d.h"

dta_dam_1d::dta_dam_1d(): numTS(0),
  year(1,1),
  month(1,1),
  day(1,1),
  Jday(1,1),
  init_year(1,1),
  init_month(1,1),
  init_day(1,1),

  InfL(1,1),
  OufL(1,1),
  OflW(1,1),
  Prec(1,1),
  Temp(1,1),
  DaiS(1,1),
  DaiG(1,1),
  OulT(1,1),
  InlT (1,1),
  EtdM(1,1),
  DamS(1,1),
  init_DamS(0.0){

  }

dta_dam_1d::~dta_dam_1d() {
  //dtor

}

dta_dam_1d::dta_dam_1d(const dta_dam_1d& other): numTS(0),
year(1,1),
month(1,1),
day(1,1),
Jday(1,1),
init_year(1,1),
init_month(1,1),
init_day(1,1),
InfL(1,1),
OufL(1,1),
OflW(1,1),
Prec(1,1),
Temp(1,1),
DaiS(1,1),
DaiG(1,1),
OulT(1,1),
InlT (1,1),
EtdM(1,1),
DamS(1,1),
init_DamS(0.0){
  numTS = other.numTS;
  year= other.year;
  month= other.month;
  day= other.day;
  Jday= other.Jday;
  init_year= other.init_year;
  init_month= other.init_month;
  init_day= other.init_day;
  InfL= other.InfL;
  OufL= other.OufL;
  OflW= other.OflW;
  Prec= other.Prec;
  Temp= other.Temp;
  DaiS= other.DaiS;
  DaiG= other.DaiG;
  EtdM= other.EtdM;
  OulT= other.OulT;
  InlT= other.InlT;
  EtdM= other.EtdM;
  DamS= other.DamS;
  init_DamS= other.init_DamS;
}

dta_dam_1d& dta_dam_1d::operator=(const dta_dam_1d& rhs) {

  if (this == &rhs) return *this;
  else {
    numTS = rhs.numTS;
    year= rhs.year;
    month= rhs.month;
    day= rhs.day;
    Jday= rhs.Jday;
    init_year= rhs.init_year;
    init_month= rhs.init_month;
    init_day= rhs.init_day;
    InfL= rhs.InfL;
    OufL= rhs.OufL;
    OflW= rhs.OflW;
    Prec= rhs.Prec;
    Temp= rhs.Temp;
    DaiS= rhs.DaiS;
    DaiG= rhs.DaiG;
    OulT= rhs.OulT;
    InlT= rhs.InlT;
    EtdM= rhs.EtdM;
    DamS= rhs.DamS;
    init_DamS= rhs.init_DamS;
  }
  return *this;
}

void dta_dam_1d::s_initStates(const hdata& initfastRes, const numberSel& init_State,const init_dStype& _Stype) {

  switch (_Stype) {
  case init_dStype::DAMS:
    init_DamS = init_State;
    break;
  }
  return ;
}

numberSel dta_dam_1d::g_initState(const init_dStype& _Stype) {

  switch (_Stype) {
  case init_dStype::DAMS:
    return init_DamS;
    break;
  }
  return 0;
}

void dta_dam_1d::s_varVal(const numberSel& dta, const unsigned& tst,const dam_ts& _tsType) {

  switch (_tsType) {
  case dam_ts::PREC:
    Prec[tst] = dta;
    break;
  case dam_ts::TEMP:
    Temp[tst] = dta;
    break;
  case dam_ts::DAMS:
    DamS[tst] = dta;
    break;
  case dam_ts::INFL:
    InfL[tst] = dta;
    break;
  case dam_ts::DAIS:
    DaiS[tst] = dta;
    break;
  case dam_ts::DAIG:
    DaiG[tst] = dta;
    break;
  case dam_ts::OUFL:
    OufL[tst] = dta;
    break;
  case dam_ts::ETDM:
    EtdM[tst] = dta;
    break;
  case dam_ts::OFLW:
    OflW[tst] = dta;
    break;
  case dam_ts::OULT:
    OulT[tst] = dta;
    break;
  case dam_ts::INLT:
    InlT[tst] = dta;
    break;
   }
  return ;


}

void dta_dam_1d::s_data(const hdata& dta,const dam_ts& _tsType, bool updateNumTS) {

  if(updateNumTS)  numTS = dta.size();

  switch (_tsType) {
  case dam_ts::PREC:
    Prec = dta;
    //    std::cout << "New precipitation --> loaded\n";
    break;
  case dam_ts::TEMP:
    Temp = dta;
    //    std::cout << "New temperature --> loaded\n";
    break;
  case dam_ts::DAMS:
    DamS = dta;
    //    std::cout << "New DamStorage --> loaded\n";
    break;
  case dam_ts::INFL:
    InfL = dta;
    //    std::cout << "New INFL --> loaded\n";
    break;
  case dam_ts::DAIS:
    DaiS = dta;
    //    std::cout << "New DAIS --> loaded\n";
    break;
  case dam_ts::DAIG:
    DaiG = dta;
    //    std::cout << "New DAIG --> loaded\n";
    break;
  case dam_ts::OUFL:
    OufL = dta;
    //    std::cout << "New OUFL --> loaded\n";
    break;
  case dam_ts::ETDM:
    EtdM = dta;
    //    std::cout << "New ETDM --> loaded\n";
    break;
  case dam_ts::OFLW:
    OflW = dta;
    //    std::cout << "New OFLW --> loaded\n";
    break;
  case dam_ts::OULT:
    OulT = dta;
    //    std::cout << "New OULT --> loaded\n";
    break;
  case dam_ts::INLT:
    InlT = dta;
    //    std::cout << "New INLT --> loaded\n";
    break;
  }
}

numberSel dta_dam_1d::g_dta(const unsigned& tst,const dam_ts& _tsType) {

  switch (_tsType) {
  case dam_ts::PREC:
    return Prec[tst];
  case dam_ts::TEMP:
    return Temp[tst];
  case dam_ts::DAMS:
    return DamS[tst];
  case dam_ts::INFL:
    return InfL[tst];
  case dam_ts::DAIS:
    return DaiS[tst];
  case dam_ts::DAIG:
    return DaiG[tst];
  case dam_ts::OUFL:
    return OufL[tst];
  case dam_ts::ETDM:
    return EtdM[tst];
  case dam_ts::OFLW:
    return OflW[tst];
  case dam_ts::OULT:
    return OulT[tst];
  case dam_ts::INLT:
    return InlT[tst];
  }

  return 0;
}

hdata dta_dam_1d::get_HbTsData(const dam_ts& _tsType) {

  switch (_tsType) {
  case dam_ts::PREC:
    return Prec;
  case dam_ts::TEMP:
    return Temp;
  case dam_ts::DAMS:
    return DamS;
  case dam_ts::INFL:
    return InfL;
  case dam_ts::DAIS:
    return DaiS;
  case dam_ts::DAIG:
    return DaiG;
  case dam_ts::OUFL:
    return OufL;
  case dam_ts::ETDM:
    return EtdM;
  case dam_ts::OFLW:
    return OflW;
  case dam_ts::OULT:
    return OulT;
  case dam_ts::INLT:
    return InlT;
  }

  hdata helpV(-9999,1);
  return helpV;
}

void dta_dam_1d::setOneTstoZero(const dam_ts& _tsType) {

  switch (_tsType) {
  case dam_ts::PREC:
    Prec = 0.0;
    break;
  case dam_ts::TEMP:
    Temp = 0.0;
    break;
  case dam_ts::DAMS:
    DamS = 0.0;
    break;
  case dam_ts::INFL:
    InfL = 0.0;
    break;
  case dam_ts::DAIS:
    DaiS = 0.0;
    break;
  case dam_ts::DAIG:
    DaiG = 0.0;
    break;
  case dam_ts::OUFL:
    OufL = 0.0;
    break;
  case dam_ts::ETDM:
    EtdM = 0.0;
    break;
  case dam_ts::OFLW:
    OflW = 0.0;
    break;
  case dam_ts::OULT:
    OulT = 0.0;
    break;
  case dam_ts::INLT:
    InlT= 0.0;
    break;
  }
  return ;
}
