#include "NuclearStructureManager.h"
#include "NMEOptions.h"
#include "Constants.h"
#include "MatrixElements.h"
#include "NuclearUtilities.h"
#include "ChargeDistributions.h"

#include "boost/algorithm/string.hpp"
#include "gsl/gsl_sf_coupling.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "NMEConfig.h"

namespace NS = NuclearStructure;
namespace NO = NS::nilsson;
namespace ME = NS::MatrixElements;
namespace CD = ChargeDistributions;

using std::cout;
using std::endl;

void ShowNMEInfo() {
  std::string author = "L. Hayen (leendert.hayen@kuleuven.be)";
  auto logger = spdlog::get("nme_results_file");
  logger->info("{:*>60}", "");
  logger->info("{:^60}", "NME v" + std::string(NME_VERSION));
  logger->info("{:^60}", "Last update: " + std::string(NME_LAST_UPDATE));
  logger->info("{:^60}", "Author: " + author);
  logger->info("{:*>60}\n", "");
}

NS::NuclearStructureManager::NuclearStructureManager() {
  InitializeLoggers();
  InitializeConstants();
}

NS::NuclearStructureManager::NuclearStructureManager(BetaType bt,
                                                     NS::Nucleus init,
                                                     NS::Nucleus fin) {
  betaType = bt;
  mother = init;
  daughter = fin;
  initialized = false;
}

void NS::NuclearStructureManager::InitializeConstants() {
  debugFileLogger->debug("Entered InitializeConstants");
  int Zd = GetNMEOpt(int, Daughter.Z);
  int Zm = GetNMEOpt(int, Mother.Z);
  int Ad = GetNMEOpt(int, Daughter.A);
  int Am = GetNMEOpt(int, Mother.A);

  if (Ad != Am) {
    consoleLogger->error(
        "ERROR: Mother and daughter mass number do not agree.");
    return;
  }
  double Rd = GetNMEOpt(double, Daughter.Radius) * 1e-15 / NATURAL_LENGTH *
              std::sqrt(5. / 3.);
  double Rm = GetNMEOpt(double, Mother.Radius) * 1e-15 / NATURAL_LENGTH *
              std::sqrt(5. / 3.);
  if (Rd == 0.0) {
    Rd = 1.2 * std::pow(Ad, 1. / 3.) * 1e-15 / NATURAL_LENGTH;
  }
  if (Rm == 0.0) {
    Rm = 1.2 * std::pow(Am, 1. / 3.) * 1e-15 / NATURAL_LENGTH;
  }
  double motherBeta2 = GetNMEOpt(double, Mother.Beta2);
  double motherBeta4 = GetNMEOpt(double, Mother.Beta4);
  double motherBeta6 = GetNMEOpt(double, Mother.Beta6);
  double daughterBeta2 = GetNMEOpt(double, Daughter.Beta2);
  double daughterBeta4 = GetNMEOpt(double, Daughter.Beta4);
  double daughterBeta6 = GetNMEOpt(double, Daughter.Beta6);
  int motherSpinParity = GetNMEOpt(int, Mother.SpinParity);
  int daughterSpinParity = GetNMEOpt(int, Daughter.SpinParity);

  double motherExcitationEn = GetNMEOpt(double, Mother.ExcitationEnergy);
  double daughterExcitationEn = GetNMEOpt(double, Daughter.ExcitationEnergy);

  std::string process = GetNMEOpt(std::string, Transition.Process);

  if (boost::iequals(process, "B+")) {
    betaType = BETA_PLUS;
  } else {
    betaType = BETA_MINUS;
  }

  if ((Zd - betaType) != Zm) {
    consoleLogger->error(
        "ERROR: Mother and daughter proton numbers cannot be coupled "
        "through process " +
        process);
    return;
  }

  SetMotherNucleus(Zm, Am, motherSpinParity, Rd, motherExcitationEn,
                   motherBeta2, motherBeta4, motherBeta6);
  SetDaughterNucleus(Zd, Ad, daughterSpinParity, Rm, daughterExcitationEn,
                     daughterBeta2, daughterBeta4, daughterBeta6);

  potential = GetNMEOpt(std::string, Computational.Potential);

  nmeResultsLogger->info("NME input overview\n{:=>30}", "");
  nmeResultsLogger->info("Using information from {}\n\n",
                         GetNMEOpt(std::string, input));
  nmeResultsLogger->info("Nuclear potential: {}", potential);
  nmeResultsLogger->info(
      "Transition from {}{} [{}/2] ({} keV) to {}{} [{}/2] ({} keV)", Am,
      utilities::atoms[int(Zm - 1)], motherSpinParity, motherExcitationEn, Ad,
      utilities::atoms[int(Zd - 1)], daughterSpinParity, daughterExcitationEn);
  debugFileLogger->debug("Leaving InitializeConstants");
}

void NS::NuclearStructureManager::InitializeLoggers() {
  std::string outputName = GetNMEOpt(std::string, output);

  /**
   * Remove result & log files if they already exist
   */
  if (std::ifstream(outputName + ".nme"))
    std::remove((outputName + ".nme").c_str());

  debugFileLogger = spdlog::get("debug_file");
  if (!debugFileLogger) {
    debugFileLogger = spdlog::basic_logger_st(
        "debug_file", GetNMEOpt(std::string, output) + ".log");
    debugFileLogger->set_level(spdlog::level::debug);
  }
  debugFileLogger->debug("Debugging logger found in NSM");
  consoleLogger = spdlog::get("console");
  if (!consoleLogger) {
    consoleLogger = spdlog::stdout_color_st("console");
    consoleLogger->set_level(spdlog::level::warn);
  }
  debugFileLogger->debug("Console logger found in NSM");
  nmeResultsLogger = spdlog::get("nme_results_file");
  if (!nmeResultsLogger) {
    nmeResultsLogger = spdlog::basic_logger_st(
        "nme_results_file", GetNMEOpt(std::string, output) + ".nme");
    nmeResultsLogger->set_level(spdlog::level::info);
    nmeResultsLogger->set_pattern("%v");
  }
  ShowNMEInfo();
  debugFileLogger->debug("NME Results logger found in NSM");
}
void NS::NuclearStructureManager::SetDaughterNucleus(int Z, int A, int dJ,
                                                     double R,
                                                     double excitationEnergy,
                                                     double beta2, double beta4,
                                                     double beta6) {
  daughter = {Z, A, dJ, R, excitationEnergy, beta2, beta4, beta6};
  initialized = false;
}

void NS::NuclearStructureManager::SetMotherNucleus(int Z, int A, int dJ,
                                                   double R,
                                                   double excitationEnergy,
                                                   double beta2, double beta4,
                                                   double beta6) {
  mother = {Z, A, dJ, R, excitationEnergy, beta2, beta4, beta6};
  initialized = false;
}

void NS::NuclearStructureManager::Initialize(std::string m, std::string p) {
  debugFileLogger->debug("Entered Initialize");
  method = m;
  // Extreme Single-particle
  if (boost::iequals(method, "ESP")) {
    potential = p;
    SingleParticleState spsi, spsf;
    int dKi, dKf;
    GetESPStates(spsi, spsf, dKi, dKf);
    for(int i = 0; i < 2; i++) {
      /* Assume one particle outside of a solid core
      */
      double robtd = std::sqrt(2*i+1.);
      AddReducedOneBodyTransitionDensity(i, robtd, dKi, dKf, spsi, spsf);
    }
    initialized = true;
  } else if (boost::iequals(method, "CUSTOMSHO")) {
    if (!NMEOptExists(Transition.DensityMatrixFile)) {
      debugFileLogger->debug(
          "Density matrix file was not specified in"
          "transition .ini file. Initializing using Method=ESP.");
      Initialize("ESP", p);
    } else {
      initialized = BuildDensityMatrixFromFile(
          GetNMEOpt(std::string, Transition.DensityMatrixFile));
    }
  }
  debugFileLogger->debug("Leaving Initialize");
}

bool NS::NuclearStructureManager::BuildDensityMatrixFromFile(
    std::string filename) {
  std::string NuShellXsuffix = ".obd";
  if (0 == filename.compare (filename.length() - NuShellXsuffix.length(), NuShellXsuffix.length(), NuShellXsuffix)){
    debugFileLogger->debug("Found NuShellX OBD file.");
    ReadNuShellXOBD(filename);
  }
  else {
    std::vector<std::vector<std::string> > dataList = GeneralUtilities::GetCSVData(filename, ",");
    for (auto const& line : dataList) {
      if (line.size() < 7) {
        consoleLogger->error("Density Matrix File {} incomplete! Aborting.",
                             filename);
        exit(EXIT_FAILURE);
      }
      int djf = atoi(line[0].c_str());
      int nf = atoi(line[1].c_str());
      int lf = atoi(line[2].c_str());
      int dji = atoi(line[3].c_str());
      int ni = atoi(line[4].c_str());
      int li = atoi(line[5].c_str());
      double obdme = atof(line[6].c_str());
      WFComp fW = {1.0, nf, lf, std::abs(djf) - 2 * lf};
      WFComp iW = {1.0, ni, li, std::abs(dji) - 2 * li};
      std::vector<WFComp> fComps = {fW};
      std::vector<WFComp> iComps = {iW};
      SingleParticleState spsf = {daughter.dJ, -1, sign(daughter.dJ), lf, nf, 0,
                                  -betaType, 0.0, fComps};
      SingleParticleState spsi = {mother.dJ, -1, sign(mother.dJ), li, ni, 0,
                                  betaType, 0.0, iComps};
      //TODO
      AddReducedOneBodyTransitionDensity(1, obdme, dji, djf, spsi, spsf);
    }
  }
  debugFileLogger->debug("Leaving BuildDensityMatrixFromFile");
  return true;
}

void NS::NuclearStructureManager::ReadNuShellXOBD(std::string filename) {
  debugFileLogger->debug("Entered ReadNuShellXOBD");
  std::string line;
  std::ifstream obdFile (filename);
  int skipHeader = 15;
  int lineNumber = 0;
  if (obdFile.is_open()) {
    while (getline(obdFile, line)) {
      lineNumber++;
      if (lineNumber > skipHeader) {
        if (line.rfind("!", 0) != 0) {
          break;
        }
      }
    }
    std::vector<std::string> stats;
    boost::algorithm::split(stats, line, boost::is_any_of(","));
    double Ji = atof(stats[0].c_str());
    double Jf = atof(stats[1].c_str());
    double Ti = atof(stats[2].c_str());
    double Tf = atof(stats[3].c_str());
    double Tip = atof(stats[4].c_str());
    double Tz = atof(stats[5].c_str());

    for(int i = 0; i < (Ji+Jf - std::abs(Ji-Jf)) + 1; i++) {
      if (getline(obdFile, line)) {
        std::vector<std::string> coupling;
        boost::algorithm::split(coupling, line, boost::is_any_of(","));
        int dJ = atoi(coupling[0].c_str());
        int ni = atoi(coupling[1].c_str());
        int nf = atoi(coupling[2].c_str());
        double ef = atof(coupling[3].c_str());
        double ei = atof(coupling[4].c_str());
        double exi = atof(coupling[5].c_str());
        double exf = atof(coupling[6].c_str());
        debugFileLogger->debug("Found coupling");
        std::vector<std::string> obd;
        while(getline(obdFile, line)) {
          if (line.rfind("   0,", 0) == 0) {
            break;
          }
          boost::algorithm::split(obd, line, boost::is_any_of(","));
          int k1 = atoi(obd[0].c_str());
          int k2 = atoi(obd[1].c_str());
          double obdMinus = std::sqrt(2.*dJ + 1.) * atof(obd[2].c_str());
          double obdPlus = std::sqrt(2.*dJ + 1.) * atof(obd[3].c_str());

          int nf = NuShellXLabels[(k1-1)*3]+1;
          int lf = NuShellXLabels[(k1-1)*3+1];
          int djf = NuShellXLabels[(k1-1)*3+2];
          int ni = NuShellXLabels[(k2-1)*3]+1;
          int li = NuShellXLabels[(k2-1)*3+1];
          int dji = NuShellXLabels[(k2-1)*3+2];

          WFComp fW = {1.0, nf, lf, std::abs(djf) - 2 * lf};
          WFComp iW = {1.0, ni, li, std::abs(dji) - 2 * li};
          std::vector<WFComp> fComps = {fW};
          std::vector<WFComp> iComps = {iW};
          SingleParticleState spsf = {djf, -1, (lf % 2 == 0) ? 1 : -1, lf, nf, 0,
                                      -betaType, 0.0, fComps};
          SingleParticleState spsi = {dji, -1, (li % 2 == 0) ? 1 : -1, li, ni, 0,
                                      betaType, 0.0, iComps};

          if (betaType == BETA_MINUS) {
            AddReducedOneBodyTransitionDensity(dJ, obdMinus, dji, djf, spsi, spsf);
            debugFileLogger->debug("Adding ROBTD with dJ = {}: {}", dJ, obdMinus);
          } else {
            AddReducedOneBodyTransitionDensity(dJ, obdPlus, dji, djf, spsi, spsf);
            debugFileLogger->debug("Adding ROBTD with dJ = {}: {}", dJ, obdPlus);
          }
        }
      }
    }
  }
  debugFileLogger->debug("Leaving ReadNuShellXOBD");
}

void NS::NuclearStructureManager::GetESPStates(SingleParticleState& spsi,
                                               SingleParticleState& spsf,
                                               int& dKi, int& dKf) {
  debugFileLogger->debug("Entered GetESPStates");
  double Vp = GetNMEOpt(double, Computational.Vproton);
  double Vn = GetNMEOpt(double, Computational.Vneutron);
  double Xn = GetNMEOpt(double, Computational.Xneutron);
  double Xp = GetNMEOpt(double, Computational.Xproton);
  double A0 = GetNMEOpt(double, Computational.SurfaceThickness);
  double VSp = GetNMEOpt(double, Computational.V0Sproton);
  double VSn = GetNMEOpt(double, Computational.V0Sneutron);

  debugFileLogger->debug("Found all Potential constants");

  int dJReqIn = mother.dJ;
  int dJReqFin = daughter.dJ;
  if (mother.A % 2 == 0) {
    dJReqIn = GetNMEOpt(int, Mother.ForcedSPSpin);
    dJReqFin = GetNMEOpt(int, Daughter.ForcedSPSpin);
  }

  double threshold = GetNMEOpt(double, Computational.EnergyMargin);

  debugFileLogger->debug("Found all spin constants");

  double V0p = Vp * (1. + Xp * (mother.A - 2. * mother.Z) / mother.A);
  double V0n = Vn * (1. - Xn * (mother.A - 2. * mother.Z) / mother.A);

  double mR = mother.R * NATURAL_LENGTH * 1e15;
  double dR = daughter.R * NATURAL_LENGTH * 1e15;

  if (boost::iequals(potential, "SHO")) {
    int ni, li, si, nf, lf, sf;
    GetESPOrbitalNumbers(ni, li, si, nf, lf, sf);
    WFComp fW = {1.0, nf, lf, sf};
    WFComp iW = {1.0, ni, li, si};
    std::vector<WFComp> fComps = {fW};
    std::vector<WFComp> iComps = {iW};
    spsf = {daughter.dJ, -1, sign(daughter.dJ), lf, nf, 0, -betaType, 0.0,
            fComps};
    spsi = {mother.dJ, -1, sign(mother.dJ), li, ni, 0, betaType, 0.0, iComps};
  } else {
    double dBeta2 = 0.0, dBeta4 = 0.0, dBeta6 = 0.0, mBeta2 = 0.0, mBeta4 = 0.0,
           mBeta6 = 0.0;
    if (boost::iequals(potential, "DWS")) {
      dBeta2 = daughter.beta2;
      dBeta4 = daughter.beta4;
      dBeta6 = daughter.beta6;
      mBeta2 = mother.beta2;
      mBeta4 = mother.beta4;
      mBeta6 = mother.beta6;
    }
    if (!GetNMEOpt(bool, Computational.ForceSpin)) {
      threshold = 0.0;
    }
    if (betaType == BETA_MINUS) {
      nmeResultsLogger->info("Proton State\n{:=>20}", "");
      spsf = NO::CalculateDeformedSPState(
          daughter.Z, 0, daughter.A, daughter.dJ, dR, dBeta2, dBeta4, dBeta6,
          V0p, A0, VSp, dJReqFin, threshold);
      nmeResultsLogger->info("Neutron State\n{:=>20}", "");
      spsi = NO::CalculateDeformedSPState(0, mother.A - mother.Z, mother.A,
                                          mother.dJ, mR, mBeta2, mBeta4, mBeta6,
                                          V0n, A0, VSn, dJReqIn, threshold);
    } else {
      nmeResultsLogger->info("Neutron State\n{:=>20}", "");
      spsf = NO::CalculateDeformedSPState(
          0, daughter.A - daughter.Z, daughter.A, daughter.dJ, dR, dBeta2,
          dBeta4, dBeta6, V0n, A0, VSn, dJReqFin, threshold);
      nmeResultsLogger->info("Proton State\n{:=>20}", "");
      spsi = NO::CalculateDeformedSPState(mother.Z, 0, mother.A, mother.dJ, mR,
                                          mBeta2, mBeta4, mBeta6, V0p, A0, VSp,
                                          dJReqIn, threshold);
    }
  }

  if (GetNMEOpt(bool, Computational.OverrideSPCoupling)) {
    dKi = std::abs(mother.dJ);
    dKf = std::abs(daughter.dJ);
  } else {
    bool reversedGhallagher = GetNMEOpt(bool, Computational.ReversedGhallagher);

    if (boost::iequals(potential, "DWS") && mother.beta2 != 0.0 &&
        daughter.beta2 != 0.0) {
      // Set Omega quantum numbers
      if (mother.A % 2 == 0) {
        int dKc = 0;
        if (((spsi.dO - spsi.lambda) == (spsf.dO - spsf.lambda) &&
             !reversedGhallagher) ||
            ((spsi.dO - spsi.lambda) != (spsf.dO - spsf.lambda) &&
             reversedGhallagher)) {
          dKc = spsi.dO + spsf.dO;
        } else {
          dKc = std::abs(spsi.dO - spsf.dO);
          if (spsi.dO > spsf.dO) {
            spsf.dO *= -1;
          } else {
            spsi.dO *= -1;
          }
        }
        if (mother.Z % 2 == 0) {
          dKi = 0;
          dKf = dKc;
        } else {
          dKf = 0;
          dKi = dKc;
        }
        if (spsf.parity * spsi.parity * dKi != mother.dJ) {
          consoleLogger->warn(
              "WARNING: Single Particle spin coupling does not match mother "
              "spin!");
          consoleLogger->warn(
              "Single Particle values. First: {}/2  Second: {}/2", spsi.dO,
              spsf.dO);
          consoleLogger->warn("Coupled spin: {} Required: {}", dKi, mother.dJ);
        }
        if (spsf.parity * spsi.parity * dKf != daughter.dJ) {
          consoleLogger->warn(
              "WARNING: Single Particle spin coupling does not match daughter "
              "spin!");
          consoleLogger->warn(
              "Single Particle values. First: {}/2  Second: {}/2", spsi.dO,
              spsf.dO);
          consoleLogger->warn("Coupled spin: {} Required: {}", dKf,
                              daughter.dJ);
        }

      } else {
        dKi = spsi.dO;
        dKf = spsf.dO;
      }
    }
  }
}

void NS::NuclearStructureManager::AddReducedOneBodyTransitionDensity(int K,
    double obdme, int dKi, int dKf, SingleParticleState spsi,
    SingleParticleState spsf) {
  ReducedOneBodyTransitionDensity robtd = {obdme, dKi, dKf, spsi, spsf};
  reducedOneBodyTransitionDensities[K].push_back(robtd);
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

double NS::NuclearStructureManager::GetESPManyParticleCoupling(
    int K, ReducedOneBodyTransitionDensity& obt) {
  double C = 0.0;

  // TODO Result from Wigner-Eckart in spin space
  int dJi = std::abs(mother.dJ);
  int dJf = std::abs(daughter.dJ);

  if (mother.A % 2 == 0) {
    int dT3i = mother.A - 2 * mother.Z;
    int dT3f = daughter.A - 2 * daughter.Z;
    int dTi = std::abs(dT3i);
    int dTf = std::abs(dT3f);
    if (NMEOptExists(Mother.Isospin)) {
      dTi = GetNMEOpt(int, Mother.Isospin);
    }
    if (NMEOptExists(NuclearPropertiesDaughterIsospin)) {
      dTf = GetNMEOpt(int, Daughter.Isospin);
    }
    /*if ((dJi + dT3i) / 2 % 2 == 0) {
      dTi = dT3i + 1;
    }
    if ((dJf + dT3f) / 2 % 2 == 0) {
      dTf = dT3f + 1;
    }*/

    debugFileLogger->debug("Isospin:");
    debugFileLogger->debug("Ti: {} Tf: {}", dTi / 2, dTf / 2);
    debugFileLogger->debug("T3: {} T3: {}", dT3i / 2, dT3f / 2);
    // Deformed transition
    if (boost::iequals(potential, "DWS") && mother.beta2 != 0.0 &&
        daughter.beta2 != 0.0) {
      if (mother.Z % 2 == 0) {
        C = 0.5 *
            std::sqrt((dJi + 1.) * (dJf + 1.) / (1. + delta(obt.dKf, 0.0))) *
            (1 + std::pow(-1., dJi / 2.)) *
            gsl_sf_coupling_3j(dJf, 2 * K, dJi, -obt.dKf, obt.dKf, 0);
      } else {
        C = 0.5 *
            std::sqrt((dJi + 1.) * (dJf + 1.) / (1. + delta(obt.dKi, 0.0))) *
            (1 + std::pow(-1., dJf / 2.)) *
            gsl_sf_coupling_3j(dJf, 2 * K, dJi, 0, -obt.dKi, obt.dKi);
      }
      // Spherical transition
    } else {
      if (mother.Z % 2 == 0) {
        // cout << "BetaType: " << betaType << endl;
        C = std::sqrt((dJi + 1.) * (dJf + 1.) * (dTi + 1.) * (dTf + 1.) /
                      (1. + delta(obt.spsi.dO, obt.spsf.dO))) *
            std::pow(-1., (dTf - dT3f) / 2.) *
            gsl_sf_coupling_3j(dTf, 2, dTi, -dT3f, -2 * betaType, dT3i) *
            gsl_sf_coupling_6j(1, dTf, (dTf + dTi) / 2, dTi, 1, 2) *
            std::sqrt(3. / 2.) * std::pow(-1., K) * 2 *
            (delta(obt.spsi.dO, obt.spsf.dO) -
             std::pow(-1., (obt.spsi.dO + obt.spsf.dO) / 2.)) *
            gsl_sf_coupling_6j(obt.spsf.dO, dJf, obt.spsi.dO, dJi, obt.spsi.dO,
                               2 * K);
      } else {
        C = std::sqrt((dJi + 1.) * (dJf + 1.) * (dTi + 1.) * (dTf + 1.) /
                      (1. + delta(obt.spsi.dO, obt.spsf.dO))) *
            std::pow(-1., (dTf - dT3f) / 2.) *
            gsl_sf_coupling_3j(dTf, 2, dTi, -dT3f, -2 * betaType, dT3i) *
            gsl_sf_coupling_6j(1, dTf, (dTf + dTi) / 2, dTi, 1, 2) *
            std::sqrt(3. / 2.) * std::pow(-1., K) * 2 *
            (1 + delta(obt.spsi.dO, obt.spsi.dO)) *
            gsl_sf_coupling_6j(obt.spsf.dO, dJf, obt.spsf.dO, dJi, obt.spsi.dO,
                               2 * K);
      }
    }
  } else {
    // cout << "Odd-A decay" << endl;
    if (boost::iequals(potential, "DWS") && mother.beta2 != 0.0 &&
        daughter.beta2 != 0.0) {
      // cout << "Deformed: " << endl;
      C = std::sqrt((dJi + 1.) * (dJf + 1.) / (1. + delta(obt.dKi, 0)) /
                    (1. + delta(obt.dKf, 0)));
    } else {
      // TODO check
      C = std::sqrt(4. * M_PI / (dJi + 1.));
      C = 1.;
    }
  }
  return C;
}

double NS::NuclearStructureManager::CalculateReducedMatrixElement(bool V, int K, int L,
                                                           int s) {
  if (!initialized) {
    Initialize(GetNMEOpt(std::string, Computational.Method),
               GetNMEOpt(std::string, Computational.Potential));
  }
  double result = 0.0;
  double nu = CD::CalcNu(mother.R * std::sqrt(3. / 5.), mother.Z);

  std::vector<ReducedOneBodyTransitionDensity> robtds = reducedOneBodyTransitionDensities[K];

  debugFileLogger->debug("Length of ROBTDS: {}", robtds.size());

  /**
  * Implementation of @f$ \langle f || \mathbf{O}_K || i \rangle = \hat{K}^{-1} \sum_{\alpha \beta} \langle \alpha || \mathbf{O}_K || \beta \rangle \langle f || [a^\dagger_\alpha \tilde{a}_\beta]_K || i \rangle @f$
  */
  for (int i = 0; i < robtds.size(); i++) {
    debugFileLogger->debug("{}", robtds[i].robtd);
    result += 1./std::sqrt(2*K+1.) * robtds[i].robtd * ME::GetReducedSingleParticleMatrixElement(
                              V, std::abs(mother.dJ) / 2., K, L, s, robtds[i].spsi,
                              robtds[i].spsf, mother.R, nu);
  }


  // /*Odd-A*/
  // int opt = 0;
  //
  // /*Even-A*/
  // if (mother.A % 2 == 0) {
  //   if (mother.Z % 2 == 1) {
  //     /*Odd-Odd*/
  //     opt = 1;
  //   } else {
  //     /*Even-Even*/
  //     opt = 2;
  //   }
  // }
  //
  // // cout << "opt: " << opt << endl;
  //
  //
  // for (int i = 0; i < oneBodyTransitions.size(); i++) {
  //   OneBodyTransition obt = oneBodyTransitions[i];
  //   if (boost::iequals(method, "ESP")) {
  //     obt.obdme = GetESPManyParticleCoupling(K, obt);
  //   }
  //   debugFileLogger->debug("OBDME: {}", obt.obdme);
  //   if (boost::iequals(potential, "DWS") && mother.beta2 != 0 &&
  //       daughter.beta2 != 0) {
  //     debugFileLogger->debug("Deformed");
  //     result += obt.obdme * ME::GetDeformedSingleParticleMatrixElement(
  //                               opt, obt.spsi, obt.spsf, V, K, L, s,
  //                               std::abs(mother.dJ), std::abs(daughter.dJ),
  //                               obt.dKi, obt.dKf, mother.R, nu);
  //   } else {
  //     result += obt.obdme * ME::GetSingleParticleMatrixElement(
  //                               V, std::abs(mother.dJ) / 2., K, L, s, obt.spsi,
  //                               obt.spsf, mother.R, nu);
  //   }
  // }
  nmeResultsLogger->info("Calculated matrix element {}M{}{}{}. Result: {}",
                         V ? "V" : "A", K, L, s, result);
  return result;
}

double NS::NuclearStructureManager::CalculateWeakMagnetism() {
  double result = 0.0;

  double gM = GetNMEOpt(double, Constants.gM);
  double gAeff = GetNMEOpt(double, Constants.gAeff);

  double VM111 = CalculateReducedMatrixElement(true, 1, 1, 1);
  double AM101 = CalculateReducedMatrixElement(false, 1, 0, 1);

  result = -std::sqrt(2. / 3.) * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV *
               mother.R / gAeff * VM111 / AM101 +
           (gM - 1) / gAeff;
  nmeResultsLogger->info("Weak magnetism final result: b/Ac = {}", result);
  return result;
}

double NS::NuclearStructureManager::CalculateInducedTensor() {
  double result = 0.0;

  double AM110 = CalculateReducedMatrixElement(false, 1, 1, 0);
  double AM101 = CalculateReducedMatrixElement(false, 1, 0, 1);

  result = 2. / std::sqrt(3.) * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV *
           mother.R * AM110 / AM101;
  nmeResultsLogger->info("Induced tensor final result: a/Ac = {}", result);
  return result;
}
