#ifndef NME_CONFIG_CONTAINER_H
#define NME_CONFIG_CONTAINER_H

#include "NHL/Containers.h"
#include "PDS/Units/GlobalSystemOfUnits.h"

#include <string>
#include <vector>

namespace BSG {

  enum class NuclearShapes {FERMI, MODGAUSS};

  struct TransitionOptions {
    int betaType = 1;
    NHL::DecayType decayType = NHL::DecayType::GAMOWTELLER;
    double mixingRatio = 0;
    double QValue;
    double atomicEnergyDeficit;
    double partialHalflife;
    double logft;
  };

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

  struct AdvancedOptions {
    bool connectSPS = false;
    NuclearShapes ESShape = NuclearShapes::FERMI;
    NuclearShapes NSShape = NuclearShapes::FERMI;
    std::vector<double> vnew, vold;
  };

  struct SpectrumCalculationOptions {
    double begin = 0.1 * keV;
    double end = 0. * keV;
    double stepSize = 0.1 * keV;
    int steps;
    bool includeNeutrino = true;
  };

  struct AllowedMatrixElements {
    double bAc;
    double dAc;
    double lambda;
  };

  struct CouplingConstants {
    double gV = 1.0;
    double gA = 1.2723;
    double gM = 4.706;
    double gP = -229.0;
  };

  struct ConfigOptions {
    CorrectionOptions correctionOptions;
    SpectrumCalculationOptions spectrumCalcOptions;
    CouplingConstants couplingConstants;
    AdvancedOptions advancedOptions;
    AllowedMatrixElements allowedME;
  };

  ConfigOptions ParseConfigFile(std::string);

  CouplingConstants ParseCouplingConstants(std::string);
  CorrectionOptions ParseCorrectionOptions(std::string);
  SpectrumCalculationOptions ParseSpectrumCalculationOptions(std::string);
  AdvancedOptions ParseAdvancedOptions(std::string);
  AllowedMatrixElements ParseAllowedMatrixElements(std::string);

  //FUTURE

  struct _FUTURE_AllowedShapeFactor {
    bool isovector;
    bool relativistic;
    bool ESDeformation;
    double bAc;
    double dAc;
    double lambda;
  };



}

#endif
