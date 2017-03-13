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

void NS::NuclearStructureManager::Initialize(std::string method, std::string pot) {
  //Extreme Single particle
  if (boost::iequals(method, "ESP")) {
    SingleParticleState spsi, spsf;
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

double NuclearStructureManager::CalculateMatrixElement(bool V, int K, int L, int s) {

}

double NuclearStructureManager::CalculateWeakMagnetism() {
  double result = 0.0;
  std::string pot = GetOpt(std::string, Computational.Potential);
  double nu = CD::CalcNu(R * std::sqrt(3. / 5.), Z);
  if (boost::iequals(pot, "SHO")) {
    //cout << "Weak magnetism SHO" << endl;
    int ni, li, si, nf, lf, sf;
    GetSPOrbitalNumbers(ni, li, si, nf, lf, sf);
    //cout << li << " " << si << " " << lf << " " << sf << endl;
    if (li == lf) {
      if (si == sf && si > 0) {
        return 1. / gA * (li + 1 + gM);
      } else if (si == sf) {
        return -1. / gA * (li - gM);
      } else {
        return 1. / 2. / gA;
      }
    }
  } else if (boost::iequals(pot, "WS")) {
    result = -std::sqrt(2. / 3.) * nucleonMasskeV / electronMasskeV * R / gA *
                 ME::GetSingleParticleMatrixElement(true, std::abs(motherSpinParity)/2., 1, 1, 1, spsi,
                                                          spsf, R, nu) /
                 ME::GetSingleParticleMatrixElement(false, std::abs(motherSpinParity)/2., 1, 0, 1, spsi,
                                                          spsf, R, nu) +
             gM / gA;
  } else if (boost::iequals(pot, "DWS")) {
    result = -std::sqrt(2. / 3.) * nucleonMasskeV / electronMasskeV * R / gA *
                 ME::CalculateDeformedSPMatrixElement(spsi, spsf, true, 1, 1, 1,
                                                      spsi.dO, spsf.dO, spsi.dK,
                                                      spsf.dK, R, nu) /
                 ME::CalculateDeformedSPMatrixElement(spsi, spsf, false, 1, 0,
                                                      1, spsi.dO, spsf.dO, spsi.dK,
                                                      spsf.dK, R, nu) +
             gM / gA;
  }
  return result;
}

double NuclearStructureManager::CalculateInducedTensor() {
  double result = 0.0;
  std::string pot = GetOpt(std::string, Computational.Potential);
  double nu = CD::CalcNu(R * std::sqrt(3. / 5.), Z);
  if (boost::iequals(pot, "SHO")) {
    int ni, li, si, nf, lf, sf;
    GetSPOrbitalNumbers(ni, li, si, nf, lf, sf);
    double dE = 2. * nu / nucleonMasskeV * (2 * (ni - nf) + li - lf);
    double rmsHO = CD::GetRMSHO(ni, li, nu);
    if (li == lf) {
      if (si == sf && si > 0) {
        return -nucleonMasskeV * dE / (2. * li + 3.) * rmsHO;
      } else if (si == sf) {
        return nucleonMasskeV * dE / (2. * li - 1.) * rmsHO;
      } else {
        return -sf * (2. * li + 1.) / 2. - dE / 2. * rmsHO;
      }
    }
  } else if (boost::iequals(pot, "WS")) {
    result =
        2. / std::sqrt(3.) * nucleonMasskeV / electronMasskeV * R *
        ME::GetSingleParticleMatrixElement(
            false, std::abs(motherSpinParity) / 2., 1, 1, 0, spsi, spsf, R, nu) / ME::GetSingleParticleMatrixElement(false, std::abs(motherSpinParity), 1, 0, 1, spsi, spsf, R, nu);
  } else if (boost::iequals(pot, "DWS")) {
    result =
        2. / std::sqrt(3.) * nucleonMasskeV / electronMasskeV * R *
        ME::CalculateDeformedSPMatrixElement(spsi, spsf, false, 1, 1, 0, spsi.dO,
                                             spsf.dO, spsi.dK, spsf.dK, R, nu) /
        ME::CalculateDeformedSPMatrixElement(spsi, spsf, false, 1, 0, 1, spsi.dO,
                                             spsf.dO, spsi.dK, spsf.dK, R, nu);
  }
  return result;
}

double NuclearStructureManager::CalculateRatioM121() {
  double M121, M101;
  std::string pot = GetOpt(std::string, Computational.Potential);
  cout << "Potential: " << pot << endl;
  double nu = CD::CalcNu(R * std::sqrt(3. / 5.), Z);
  if (boost::iequals(pot, "SHO")) {
    int ni, li, si, nf, lf, sf;
    GetSPOrbitalNumbers(ni, li, si, nf, lf, sf);
    M121 = ME::GetSingleParticleMatrixElement(
        false, std::abs(motherSpinParity) / 2., 1, 2, 1, ni, nf, li, lf, si, sf,
        R, nu);
    M101 = ME::GetSingleParticleMatrixElement(
        false, std::abs(motherSpinParity) / 2., 1, 0, 1, ni, nf, li, lf, si, sf,
        R, nu);
  } else if (boost::iequals(pot, "WS")) {
    M121 = ME::GetSingleParticleMatrixElement(
        false, std::abs(motherSpinParity) / 2., 1, 2, 1, spsi, spsf, R, nu);
    M101 = ME::GetSingleParticleMatrixElement(
        false, std::abs(motherSpinParity) / 2., 1, 0, 1, spsi, spsf, R, nu);
  } else if (boost::iequals(pot, "DWS")) {
    M121 = ME::CalculateDeformedSPMatrixElement(
        spsi, spsf, false, 1, 2, 1, spsi.dO, spsf.dO, spsi.dK, spsf.dK, R, nu);
    M101 = ME::CalculateDeformedSPMatrixElement(
        spsi, spsf, false, 1, 0, 1, spsi.dO, spsf.dO, spsi.dK, spsf.dK, R, nu);
  }
  cout << "M121: " << M121 << " M101: " << M101 << endl;
  fc1 = gA * M101;

  return M121 / M101;
}
