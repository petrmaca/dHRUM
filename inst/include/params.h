#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <valarray>
#include <vector>
#include <utility> // pair
#include <string> //string
#include <list>

#include "numberSel.h"
#include "parStructSels.h"


class params {
 public:
  params();
  virtual ~params();
  params(const params& other);
  params& operator=(const params& other);

  void s_params(const numberSel& dta,par_HRUtype _parType);//!< The setting of model parameters
  numberSel g_par(const par_HRUtype& _parType);//!< The getting of model parameters
  numberSel g_par_up(const par_HRUtype& _parType);//!< Getting upper bounds on model parameters
  numberSel g_par_low(const par_HRUtype& _parType);//!< Getting lower bounds on model parameters

  unsigned g_numFastRes();//!< Get the number of fast reservoirs
  void s_numFastRes(const unsigned& numRes);//!< set number of fast reservoirs


  void s_parLoadToCalib(std::vector<std::pair<numberSel,par_HRUtype>>& parsToLoad);
  void s_params(const std::pair <numberSel, par_HRUtype>& parDta);
  void s_default();

  void p_param();
  void PDM_boundary_update(); //!< Adjusts the upper and lower parameters cmin and cmax when using the PDM model

  unsigned g_numPars();//!< Get the number of parameters
  void current_param(gs_STORtype gs_STORAGE,soil_STORtype soil_STORAGE,interception_STORtype intrc_STORAGE,surface_STORtype srfs_STORAGE,fast_Response fast_RESP);
  std::vector<std::string> par_HRUtype_to_string(std::list<par_HRUtype> par_list);
  void print_par_list(std::list<par_HRUtype> par_list);

  std::list<par_HRUtype> Current_parameter_list;
  // std::vector<std::string> Current_parameter_string;
  std::vector<double> Current_parameter_val;
  std::vector<double> Current_upparameter_val;
  std::vector<double> Current_lowparameter_val;

  unsigned g_sizeVecNamesPars();


 protected:

 private:
  //!< The parameters for single HMunit
  unsigned numPars;//!<  The number of model parameters
  pdata pars;//!< The values of parameters on given HRU
  pdata up_pars;//!< The upper bounds of parameters
  pdata low_pars;//!< The lower bound of parameters
  unsigned numFastRes;//!< The number of fast runoff reservoirs, states are implemented in data_HB_1d class

  std::vector<std::string> Current_parameter_string;//!< The vector of param names
/*
  //srfs_STORAGE
  std::list<par_HRUtype> L_SurfaceAll = { par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::TETR, par_HRUtype::RETCAP };
  std::list<par_HRUtype> L_SurfacePRTL = {  };
  std::list<par_HRUtype> L_Wetland = {  };
  //intrc_STORAGE
  std::list<par_HRUtype> L_Rutter_Gash = { par_HRUtype::CDIV, par_HRUtype::SDIV, par_HRUtype::CAN_ST, par_HRUtype::STEM_ST,par_HRUtype::CSDIV };
  //gs_STORAGE
  std::list<par_HRUtype> L_LIN_RES = { par_HRUtype::KS, par_HRUtype::ADIV };
  std::list<par_HRUtype> L_LINL_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::L };
  std::list<par_HRUtype> L_LINBY_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::D_BYPASS };
  std::list<par_HRUtype> L_POW_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::B_EXP };
  std::list<par_HRUtype> L_EXP_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::B_EXP };
  std::list<par_HRUtype> L_LIN_2SE = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::KS2 };
  std::list<par_HRUtype> L_LIN_2PA = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::KS2, par_HRUtype::ALPHA };
  std::list<par_HRUtype> L_FLEX_RES = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::KS2, par_HRUtype::THR };
  std::list<par_HRUtype> L_EXP_LOG = { par_HRUtype::KS, par_HRUtype::ADIV, par_HRUtype::B_EXP };
  //soil_STORAGE
  std::list<par_HRUtype> L_PDM = { par_HRUtype::B_SOIL, par_HRUtype::C_MAX, par_HRUtype::B_EVAP, par_HRUtype::SMAXpdm, par_HRUtype::CMIN };
  std::list<par_HRUtype> L_COLLIE_V2 = { par_HRUtype::KF, par_HRUtype::FC, par_HRUtype::FOREST_FRACT, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_NEW_ZEALAND = { par_HRUtype::KF, par_HRUtype::FC, par_HRUtype::FOREST_FRACT, par_HRUtype::KF2, par_HRUtype::KF_NONLIN, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_GR4J = { par_HRUtype::SMAX};
  std::list<par_HRUtype> L_SBROOK_V1 = { par_HRUtype::KF, par_HRUtype::FC, par_HRUtype::FOREST_FRACT, par_HRUtype::KF_NONLIN, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_HILLSLOPE = { par_HRUtype::KF, par_HRUtype::KF_NONLIN, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_PLATEAU = { par_HRUtype::C, par_HRUtype::INFR_MAX, par_HRUtype::RF, par_HRUtype::WP, par_HRUtype::SMAX};
  std::list<par_HRUtype> L_PDM2 = { par_HRUtype::B_SOIL, par_HRUtype::C_MAX };
  //fast_response
  std::list<par_HRUtype> L_SerialCascadeLinRes = { par_HRUtype::KFR, par_HRUtype::ADIV };
  std::list<par_HRUtype> L_SerialLinResGWGros = { par_HRUtype::KFR, par_HRUtype::ADIV,par_HRUtype::RBEI };
  std::list<par_HRUtype> L_SerialLinResSoilSois = { par_HRUtype::KFR, par_HRUtype::ADIV,par_HRUtype::RBAI };
  std::list<par_HRUtype> L_SerialLinResGWGrosSoilSois = { par_HRUtype::KFR, par_HRUtype::ADIV,par_HRUtype::RBEI, par_HRUtype::RBAI };
  //other
  std::list<par_HRUtype> L_interception_snow = { par_HRUtype::TETR };
  std::list<par_HRUtype> L_snow_melt = { par_HRUtype::DDFA, par_HRUtype::TMEL };
*/
};

#endif // PARAMS_H

/*
#include <iostream>
#include <list>
#include <vector>
#include <string>

// Define the core enum class for the traffic signals
enum class TrafficLight {
  Red,
  Yellow,
  Green,
  FlashingYellow,
  Off
};

// Define an enum to represent the configuration mode of the intersection
enum class ControllerMode {
  Standard,
  PedestrianCrossing,
  NightMode
};

// Helper function to convert enum to string for clean output
std::string status_to_string(TrafficLight status) {
  switch (status) {
  case TrafficLight::Red: return "Red";
  case TrafficLight::Yellow: return "Yellow";
  case TrafficLight::Green: return "Green";
  case TrafficLight::FlashingYellow: return "Flashing Yellow";
  case TrafficLight::Off: return "Off";
  }
  return "Unknown";
}

class IntersectionController {
private:
  // **NEW CONST MEMBER:** This list is constant and holds the default fixed sequence.
  // It MUST be initialized in the Member Initializer List.
  const std::list<TrafficLight> standard_sequence;

  // Private member: a list containing TrafficLight enum types (mutable, based on mode)
  std::list<TrafficLight> signal_sequence;
  ControllerMode mode;

  // Helper function to generate the correct sequence based on the mode
  // This is called by the initializer list of the parameterized constructor.
  static std::list<TrafficLight> create_sequence(ControllerMode m) {
    switch (m) {
    case ControllerMode::Standard:
      return {
      TrafficLight::Red,
      TrafficLight::Yellow,
      TrafficLight::Green,
      TrafficLight::Yellow
    };
    case ControllerMode::PedestrianCrossing:
      // A mode that includes a longer Red phase for pedestrians
      return {
      TrafficLight::Red,
      TrafficLight::Off, // Brief pause
      TrafficLight::Red
    };
    case ControllerMode::NightMode:
      // A mode that only uses flashing yellow for caution
      return {
      TrafficLight::FlashingYellow,
      TrafficLight::FlashingYellow
    };
    }
    return {}; // Default/empty sequence
  }

  // Helper to create the fixed standard sequence for the const member
  static std::list<TrafficLight> create_standard_sequence() {
    return {
    TrafficLight::Red,
    TrafficLight::Yellow,
    TrafficLight::Green,
    TrafficLight::Yellow
  };
  }

public:
  // 1. Default Constructor (Uses the previous fixed sequence)
  IntersectionController() :
  // Initializing the new const member first
  standard_sequence(create_standard_sequence()),
  // Initializing the mutable members after
  mode(ControllerMode::Standard),
  signal_sequence(standard_sequence) // Can initialize the mutable list from the const list
  {
    std::cout << "Controller initialized in STANDARD mode with " << signal_sequence.size() << " signals." << std::endl;
    std::cout << "Const list size: " << standard_sequence.size() << std::endl;
  }

  // 2. Parameterized Constructor (Initializes subset based on mode)
  IntersectionController(ControllerMode initial_mode) :
  // Initializing the new const member first (must be done!)
  standard_sequence(create_standard_sequence()),
  // Initializing the mutable members after
  mode(initial_mode),
  // Crucial step: The initializer list calls a static helper function
  // to return the appropriate list subset based on the constructor argument.
  signal_sequence(create_sequence(initial_mode))
  {
    std::cout << "Controller initialized in CUSTOM mode (Mode " << static_cast<int>(initial_mode) << ") with " << signal_sequence.size() << " signals." << std::endl;
    std::cout << "Const list size: " << standard_sequence.size() << std::endl;
  }

  // Public method to safely display the sequence
  void display_sequence() const {
    std::cout << "\n[Mode: " << (mode == ControllerMode::Standard ? "Standard" : (mode == ControllerMode::PedestrianCrossing ? "Pedestrian" : "Night")) << "] Sequence: ";
    bool first = true;
    for (const auto& status : signal_sequence) {
      if (!first) {
        std::cout << " -> ";
      }
      std::cout << status_to_string(status);
      first = false;
    }
    std::cout << "\n" << std::endl;
  }

  // Read-only getter for the constant list
  const std::list<TrafficLight>& get_standard_sequence() const {
    return standard_sequence;
  }
};

int main() {
  // 1. Initialize using the Default Constructor (Standard Mode)
  IntersectionController standard_controller;
  standard_controller.display_sequence();

  // Demonstrate accessing the const member (read-only)
  std::cout << "Checking const sequence via getter. First item is: " << status_to_string(standard_controller.get_standard_sequence().front()) << std::endl;

  // 2. Initialize a Pedestrian Crossing controller instance
  IntersectionController pedestrian_controller(ControllerMode::PedestrianCrossing);
  pedestrian_controller.display_sequence();

  // 3. Initialize a Night Mode controller instance
  IntersectionController night_controller(ControllerMode::NightMode);
  night_controller.display_sequence();

  return 0;
}
*/
