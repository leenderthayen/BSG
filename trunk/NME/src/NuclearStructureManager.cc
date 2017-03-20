#include "NuclearStructureManager.h"
#include "OptionContainer.h"
#include "Constants.h"
#include "MatrixElements.h"
#include "NuclearUtilities.h"
#include "ChargeDistributions.h"

#include "boost/algorithm/string.hpp"
#include "gsl/gsl_sf_coupling.h"

namespace NS = NuclearStructure;
namespace NO = NS::nilsson;
namespace ME = NS::MatrixElements;
namespace CD = ChargeDistributions;

NS::NuclearStructureManager::NuclearStructureManager(BetaType bt,
                                                     NS::Nucleus init,
                                                     NS::Nucleus fin) {
  betaType = bt;
  mother = init;
  daughter = fin;
}

void NS::NuclearStructureManager::SetDaughterNucleus(int Z, int A, int dJ,
                                                     double R,
                                                     double excitationEnergy,
                                                     double beta2,
                                                     double beta4) {
  daughter = {Z, A, dJ, R, excitationEnergy, beta2, beta4};
}

void NS::NuclearStructureManager::SetMotherNucleus(int Z, int A, int dJ,
                                                   double R,
                                                   double excitationEnergy,
                                                   double beta2, double beta4) {
  mother = {Z, A, dJ, R, excitationEnergy, beta2, beta4};
}

void NS::NuclearStructureManager::Initialize(std::string m,
                                             std::string p) {
  method = m;
  potential = p;
  // Extreme Single particle
  if (boost::iequals(method, "ESP")) {
    SingleParticleState spsi, spsf;
    int dKi, dKf;
    double obdme;
    GetESPStates(spsi, spsf, potential, dKi, dKf, obdme);
    AddOneBodyTransition(obdme, dKi, dKf, spsi, spsf);
  }
}

void NS::NuclearStructureManager::GetESPStates(SingleParticleState& spsi,
                                               SingleParticleState& spsf,
                                               std::string potential, int& dKi,
                                               int& dKf, double& obdme) {
  double V0p = GetOpt(double, Computational.V0proton);
  double V0n = GetOpt(double, Computational.V0neutron);
  double A0 = GetOpt(double, Computational.SurfaceThickness);
  double VSp = GetOpt(double, Computational.V0Sproton);
  double VSn = GetOpt(double, Computational.V0Sneutron);

  double mR = mother.R * NATLENGTH * 1e15;
  double dR = daughter.R * NATLENGTH * 1e15;

  if (boost::iequals(potential, "SHO")) {
    int ni, li, si, nf, lf, sf;
    GetESPOrbitalNumbers(ni, li, si, nf, lf, sf);
    WFComp fW = {1.0, nf, lf, sf};
    WFComp iW = {1.0, ni, li, si};
    std::vector<WFComp> fComps = {fW};
    std::vector<WFComp> iComps = {iW};
    spsf = {daughter.dJ, -1, sign(daughter.dJ), -betaType, 0.0, fComps};
    spsi = {mother.dJ, -1, sign(mother.dJ), betaType, 0.0, iComps};
  } else {
    double dBeta2 = 0.0, dBeta4 = 0.0, mBeta2 = 0.0, mBeta4 = 0.0;
    if (boost::iequals(potential, "DWS")) {
      dBeta2 = daughter.beta2;
      dBeta4 = daughter.beta4;
      mBeta2 = mother.beta2;
      mBeta4 = mother.beta4;
    }
    if (betaType == BETA_MINUS) {
      spsf =
          NO::CalculateDeformedSPState(daughter.Z, 0, daughter.A, daughter.dJ,
                                       dR, dBeta2, dBeta4, V0p, A0, VSp);
      spsi = NO::CalculateDeformedSPState(0, mother.A - mother.Z, mother.A,
                                          mother.dJ, mR, mBeta2, mBeta4, V0n,
                                          A0, VSn);
    } else {
      spsf = NO::CalculateDeformedSPState(0, daughter.A - daughter.Z,
                                          daughter.A, daughter.dJ, dR, dBeta2,
                                          dBeta4, V0n, A0, VSn);
      spsi = NO::CalculateDeformedSPState(mother.Z, 0, mother.A, mother.dJ, mR,
                                          mBeta2, mBeta4, V0p, A0, VSp);
    }
  }

  // Set Omega quantum numbers
  if (mother.A % 2 == 0) {
    int dTi, dTf;
    int dT3i = std::abs(mother.A - 2 * mother.Z);
    int dT3f = std::abs(daughter.A - 2 * daughter.Z);
    if ((mother.dJ + dT3i) / 2 % 2 == 0) {
      dTi = dT3i + 1;
    }
    if ((daughter.dJ + dT3f) / 2 % 2 == 0) {
      dTf = dT3f + 1;
    }
    // TODO Implement Gallaghar coupling rules
    // Even-Even ---> Odd-Odd transition
    if (mother.Z % 2 == 0) {
      dKi = 0;
      dKf = spsi.dO + spsf.dO;
      // Deformed transition
      if (boost::iequals(potential, "DWS") && mother.beta2 != 0.0 &&
          daughter.beta2 != 0.0) {
        obdme = 0.5 * std::sqrt((mother.dJ + 1.) * (daughter.dJ + 1.) / 2. /
                                (1. + delta(dKf, 0.0))) *
                (1 + std::pow(-1., (dKf + spsi.dO + spsf.dO) / 2.));
        // Spherical transition
      } else {
        obdme =
            std::sqrt((mother.dJ + 1.) * (daughter.dJ + 1.) * (dTi + 1.) *
                      (dTf + 1.) / 2. /
                      (1. + delta(spsi.dO, spsf.dO))) *
            std::pow(-1., (dTf - dT3f) / 2.) *
            gsl_sf_coupling_3j(dTf, 2, dTi, dT3f, 2 * betaType, dT3i) *
            gsl_sf_coupling_6j(1, dTf, 1, dTi, 1, 2) * std::sqrt(3. / 2.) * 2 *
            (delta(spsi.dO, spsf.dO) -
             std::pow(-1., (spsi.dO + spsf.dO) / 2.));
      }
      // Odd-Odd ---> Even-Even transition
    } else {
      dKf = 0;
      dKi = spsi.dO + spsf.dO;
      if (boost::iequals(potential, "DWS") && mother.beta2 != 0.0 &&
          daughter.beta2 != 0.0) {
        obdme = 0.5 * std::sqrt((mother.dJ + 1.) * (daughter.dJ + 1.) / 2. /
                                (1. + delta(dKi, 0.0))) *
                (1 + std::pow(-1., daughter.dJ / 2.));
      } else {
        obdme =
            std::sqrt((mother.dJ + 1.) * (daughter.dJ + 1.) * (dTi + 1.) *
                      (dTf + 1.) / 2. /
                      (1. + delta(spsi.dO, spsf.dO))) *
            std::pow(-1., (dTf - dT3f) / 2.) *
            gsl_sf_coupling_3j(dTf, 2, dTi, dT3f, 2 * betaType, dT3i) *
            gsl_sf_coupling_6j(1, dTf, 1, dTi, 1, 2) * std::sqrt(3. / 2.) * 2 *
            (1 + delta(spsi.dO, spsi.dO));
      }
    }
  } else {
    dKi = spsi.dO;
    dKf = spsf.dO;
    if (boost::iequals(potential, "DWS") && mother.beta2 != 0.0 &&
        daughter.beta2 != 0.0) {
      obdme = std::sqrt((mother.dJ + 1.) * (daughter.dJ + 1.) /
                        (1. + delta(dKi, 0)) / (1. + delta(dKf, 0)));
    } else {
      // TODO check
      obdme = std::sqrt(4. * M_PI / (mother.dJ + 1.));
    }
  }
}

void NS::NuclearStructureManager::GetESPOrbitalNumbers(int& ni, int& li,
                                                       int& si, int& nf,
                                                       int& lf, int& sf) {
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

double NS::NuclearStructureManager::CalculateMatrixElement(bool V, int K, int L,
                                                           int s) {
  double result = 0.0;
  double nu = CD::CalcNu(mother.R * std::sqrt(3. / 5.), mother.Z);
  int opt = 0;

  if (mother.A % 2 == 0) {
    if (mother.Z % 2 == 1) {
      opt = 1;
    } else {
      opt = 2;
    }
  }

  for (int i = 0; i < oneBodyTransitions.size(); i++) {
    OneBodyTransition obt = oneBodyTransitions[i];
    if (boost::iequals(potential, "DWS") && mother.beta2 != 0 &&
        daughter.beta2 != 0) {
      result += obt.obdme * ME::GetDeformedSingleParticleMatrixElement(opt,
                                obt.spsi, obt.spsf, V, K, L, s, mother.dJ,
                                daughter.dJ, obt.dKi, obt.dKf, mother.R, nu);
    } else {
      result += obt.obdme * ME::GetSingleParticleMatrixElement(
                                V, std::abs(mother.dJ) / 2., K, L, s, obt.spsi,
                                obt.spsf, mother.R, nu);
    }
  }
  return result;
}

double NS::NuclearStructureManager::CalculateWeakMagnetism() {
  double result = 0.0;

  double VM111 = CalculateMatrixElement(true, 1, 1, 1);
  double AM101 = CalculateMatrixElement(false, 1, 0, 1);

  result = -std::sqrt(2. / 3.) * nucleonMasskeV / electronMasskeV * mother.R /
               gA * VM111 / AM101 +
           gM / gA;
  return result;
}

double NS::NuclearStructureManager::CalculateInducedTensor() {
  double result = 0.0;

  double AM110 = CalculateMatrixElement(false, 1, 1, 0);
  double AM101 = CalculateMatrixElement(false, 1, 0, 1);

  result = 2. / std::sqrt(3.) * nucleonMasskeV / electronMasskeV * mother.R *
           AM110 / AM101;
  return result;
}
