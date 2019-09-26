#ifndef NME_CONFIG_CONTAINER_H
#define NME_CONFIG_CONTAINER_H

#include "PDS/Units/GlobalSystemOfUnits.h"

#include <string>

namespace BSG {

  struct CorrectionOptions{
    bool phaseSpace = true;
    bool FermiFunction = true;
    bool ESFiniteSize = true;
    bool shapeFactor = true;
    bool isovector = true;
    bool relativistic = true;
    bool ESDeformation = true;
    bool ESFermi = true;
    bool CoulombRecoil = true;
    bool radiative = true;
    bool kinRecoil = true;
    bool atomicScreen = true;
    bool atomicExchange = true;
    bool atomicMismatch = true;
  };

  struct SpectrumCalculationOptions {
    double begin = 0.1 * keV;
    double end = 0. * keV;
    double stepSize = 0.1 * keV;
    int steps;
    bool includeNeutrino = true;
  };

  struct CouplingConstants {
    double gV = 1.0;
    double gA = 1.2723;
    double gM = 4.706;
    double gP = -229.0;
  };

  struct ConfigOptions {
    ROBTDMethod robtdMethod;
    ESPOptions espOptions;
    WSPotentialOptions wsPotentialOptions;
    CouplingConstants couplingConstants;
  };

  ConfigOptions ParseConfigFile(std::string);

  CouplingConstants ParseCouplingConstants(std::string);
  WSPotentialOptions ParseWSPotentialOptions(std::string);
  ESPOptions ParseESPOptions(std::string);

}

#endif
