#include "data_dam_1d.h"

data_dam_1d::data_dam_1d(): numTS(0),
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
  Sois(1,1),
  Gros(1,1),
  OulT(1,1),
  InlT (1,1),
  EtdM(1,1),
  DamS(1,1),
  damMax(0.0),
  MRF(0.0),
  damArea(0.0),
  damLeak(0.0),
  init_DamS(0.0){

  }

data_dam_1d::~data_dam_1d() {
  //dtor

}

data_dam_1d::data_dam_1d(const data_dam_1d& other): numTS(0),
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
Sois(1,1),
Gros(1,1),
EtdM(1,1),
OulT(1,1),
InlT (1,1),
DamS(1,1),
damMax(0.0),
MRF(0.0),
damArea(0.0),
damLeak(0.0),
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
  Sois= other.Sois;
  Gros= other.Gros;
  EtdM= other.EtdM;
  OulT= other.OulT;
  InlT= other.InlT;
  EtdM= other.EtdM;
  DamS= other.DamS;
  damMax= other.damMax;
  MRF= other.MRF;
  damArea= other.damArea;
  damLeak= other.damLeak;
  init_DamS= other.init_DamS;
}

data_dam_1d& data_dam_1d::operator=(const data_dam_1d& rhs) {

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
    Sois= rhs.Sois;
    Gros= rhs.Gros;
    EtdM= rhs.EtdM;
    OulT= rhs.OulT;
    InlT= rhs.InlT;
    EtdM= rhs.EtdM;
    DamS= rhs.DamS;
    damMax= rhs.damMax;
    MRF= rhs.MRF;
    damArea= rhs.damArea;
    damLeak= rhs.damLeak;
    init_DamS= rhs.init_DamS;
  }
  return *this;
}
