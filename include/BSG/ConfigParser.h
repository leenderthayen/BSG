#ifndef NME_CONFIG_CONTAINER_H
#define NME_CONFIG_CONTAINER_H

#include "NHL/Containers.h"
#include "PDS/Units/GlobalSystemOfUnits.h"
#include "PDS/Core/ParticleDefinition.h"

#include "CLI11.hpp"

#include <string>
#include <vector>

namespace BSG {

  struct TransitionOptions {
    PDS::core::Particle initNucleus, finalNucleus;
    NHL::BetaType betaType = NHL::BETA_MINUS;
    NHL::BetaDecayType decayType = NHL::BetaDecayType::GAMOWTELLER;
    double mixingRatio = 0;
    double QValue = 0.0 * keV;
    double atomicEnergyDeficit = 0.0 * keV;
    double partialHalflife = 0.0 * s;
    double logft = 0.;
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
    bool atomicScreening = true;
    bool atomicExchange = false;
    bool atomicMismatch = true;
  };

  struct AdvancedOptions {
    bool connectSPS = false;
    NHL::NuclearShapes ESShape = NHL::NuclearShapes::FERMI;
    NHL::NuclearShapes NSShape = NHL::NuclearShapes::FERMI;
    std::vector<double> vnew, vold;
    double modGaussFit = 5;
  };

  struct SpectrumCalculationOptions {
    double begin = 0.1 * keV;
    double end = 0. * keV;
    double stepSize = 10 * keV;
    int steps = 0;
    bool includeNeutrino = true;
  };

  struct AllowedMatrixElements {
    double bAc = 0.;
    double dAc = 0.;
    double lambda = 0.;
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

  TransitionOptions ParseTransitionOptions(std::string, int argc = 0, const char** argv = nullptr);
  void SetTransitionOptions(CLI::App&, TransitionOptions&);

  ConfigOptions ParseConfigOptions(std::string, int argc = 0, const char** argv = nullptr);

  void SetCouplingOptions(CLI::App&, CouplingConstants&);
  void SetCorrectionOptions(CLI::App&, CorrectionOptions&);
  void SetSpectrumCalculationOptions(CLI::App&, SpectrumCalculationOptions&);
  void SetAdvancedOptions(CLI::App&, AdvancedOptions&);
  void SetAllowedMatrixElementOptions(CLI::App&, AllowedMatrixElements&);

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
