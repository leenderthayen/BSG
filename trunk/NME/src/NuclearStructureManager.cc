#include "NuclearStructureManager"
#include "OptionContainer.h"
#include "Constants.h"

#include "boost/algorithm/string.hpp"

namespace NS = NuclearStructure;
namespace NO = NuclearStrucrure::nilsson;

NS::NuclearStructureManager::NuclearStructureManager(BetaType bt, NS::Nucleus init, NS::Nucleus fin) {
  betaType = bt;
  mother = init;
  daughter = fin;
}

void NS::NuclearStructureManager::SetDaughterNucleus(int Z, int A, int dJ, double R, double excitationEnergy, double beta2, double beta4) {
  daughter = {Z, A, dJ, R, excitationEnergy, beta2, beta4};
}

void NS::NuclearStructureManager::SetMotherNucleus(int Z, int A, int dJ, double R, double excitationEnergy, double beta2, double beta4) {
  mother = {Z, A, dJ, R, excitationEnergy, beta2, beta4};
}

void NS::NuclearStructureManager::Initialize(std::string method) {
  //Extreme Single particle
  if (boost::iequals(method, "ESP")) {
    SingleParticleState spsi, spsf;
    std::string pot = GetOpt(std::string, Computational.Potential);
    GetESPStates(spsi, spsf, pot);
    
  }
}

void NS::NuclearStructureManager::GetESPStates(SingleParticleState& spsi, SingleParticleState& spsf, std::string pot) {
  double V0p = GetOpt(double, Computational.V0proton);
  double V0n = GetOpt(double, Computational.V0neutron);
  double A0 = GetOpt(double, Computational.SurfaceThickness);
  double VSp = GetOpt(double, Computational.V0Sproton);
  double VSn = GetOpt(double, Computational.V0Sneutron);

  double mR = mother.R * NATLENGTH * 1e15;
  double dR = daughter.R * NATLENGTH * 1e15;

  if (boost::iequals(pot, "SHO")) {
    int ni, li, si, nf, lf, sf;
    GetESPOrbitalNumbers(ni, li, si, nf, lf, sf);
    WFComp fW = {1.0, nf, lf, sf};
    WFComp iW = {1.0, ni, li, si};
    spsf = {}
  } else if (boost::iequals(pot, "WS")) {
    if (betaType == BETA_MINUS) {
      spsf =
          NO::CalculateDeformedState(daughter.Z, 0, daughter.A, daughter.dJ, dR, 0.0, 0.0, V0p, A0, VSp);
      spsi = nilsson::CalculateDeformedState(0, mother.A - mother.Z, mother.A, mother.dJ, mR, 0.0,
                                             0.0, V0n, A0, VSn);
    } else {
      spsf = nilsson::CalculateDeformedState(0, A - Z, A, daughterSpinParity, r, 0.0, 0.0, V0n, A0,
                                             VSn);
      spsi = nilsson::CalculateDeformedState(Z - fBetaType, 0, A, motherSpinParity, r, 0.0, 0.0,
                                             V0p, A0, VSp);
    }
  } else if (boost::iequals(pot, "DWS")) {
    if (fBetaType == BETA_MINUS) {
      spsf = nilsson::CalculateDeformedState(Z, 0, A, daughterSpinParity, r, daughterBeta2,
                                             daughterBeta4, V0p, A0, VSp);
      spsi = nilsson::CalculateDeformedState(0, A - (Z - fBetaType), A, motherSpinParity, r,
                                             motherBeta2, motherBeta4, V0n, A0,
                                             VSn);
    } else {
      cout << "Final state" << endl;
      spsf = nilsson::CalculateDeformedState(0, A - Z, A, daughterSpinParity, r, daughterBeta2,
                                             daughterBeta4, V0n, A0, VSn);
      cout << "Initial state: " << endl;
      spsi = nilsson::CalculateDeformedState(
          Z - fBetaType, 0, A, motherSpinParity, r, motherBeta2, motherBeta4, V0p, A0, VSp);
    }
  }
}

void NS::NuclearStructureManager::GetESPOrbitalNumbers(int& ni, int& li, int& si, int& nf, int& lf,
                                    int& sf) {
  std::vector<int> occNumbersInit, occNumbersFinal;
  if (betaType == BETA_MINUS) {
    occNumbersInit = utilities::GetOccupationNumbers(mother.A - mother.Z);
    occNumbersFinal = utilities::GetOccupationNumbers(daughter.Z);
  } else {
    occNumbersInit = utilities::GetOccupationNumbers(mother.Z);
    occNumbersFinal = utilities::GetOccupationNumbers(daughter.A - daughter.Z);
  }
  ni = occNumbersInit[occNumbersInit.size() - 1 - 3];
  li = occNumbersInit[occNumbersInit.size() - 1 - 2];
  si = occNumbersInit[occNumbersInit.size() - 1 - 1];
  nf = occNumbersFinal[occNumbersFinal.size() - 1 - 3];
  lf = occNumbersFinal[occNumbersFinal.size() - 1 - 2];
  sf = occNumbersFinal[occNumbersFinal.size() - 1 - 1];
}
