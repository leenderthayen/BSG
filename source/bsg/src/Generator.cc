#include "Generator.h"

#include "OptionContainer.h"
#include "ChargeDistributions.h"
#include "Constants.h"
#include "Utilities.h"
#include "SpectralFunctions.h"

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <chrono>

#include "boost/algorithm/string.hpp"

#include "BSGConfig.h"

namespace SF = SpectralFunctions;
namespace CD = ChargeDistributions;
namespace NS = NuclearStructure;

using std::cout;
using std::endl;

void ShowBSGInfo() {
  std::string author = "L. Hayen (leendert.hayen@kuleuven.be)";
  auto logger = spdlog::get("BSG_results_file");
  logger->info("{:*>60}", "");
  logger->info("{:^60}", "BSG v" + std::string(BSG_VERSION));
  logger->info("{:^60}", "Last update: " + std::string(BSG_LAST_UPDATE));
  logger->info("{:^60}", "Author: " + author);
  logger->info("{:*>60}\n", "");
}

Generator::Generator() {
  auto start = std::chrono::steady_clock::now();
  InitializeLoggers();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "Loggers milliseconds since start: " << elapsed.count() << "\n";
  InitializeConstants();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "Constants milliseconds since start: " << elapsed.count() << "\n";
  InitializeShapeParameters();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "Shape milliseconds since start: " << elapsed.count() << "\n";
  InitializeL0Constants();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "L0 milliseconds since start: " << elapsed.count() << "\n";
  LoadExchangeParameters();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "Exchange milliseconds since start: " << elapsed.count() << "\n";

  InitializeNSMInfo();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "NSM milliseconds since start: " << elapsed.count() << "\n";
}

Generator::~Generator() { delete nsm; }

void Generator::InitializeLoggers() {
  std::string outputName = GetBSGOpt(std::string, output);

  /**
   * Remove result & log files if they already exist
   */
  if (std::ifstream(outputName + ".log")) std::remove((outputName + ".log").c_str());
  if (std::ifstream(outputName + ".raw")) std::remove((outputName + ".raw").c_str());
  if (std::ifstream(outputName + ".txt")) std::remove((outputName + ".txt").c_str());

  debugFileLogger = spdlog::get("debug_file");
  if (!debugFileLogger) {
    debugFileLogger = spdlog::basic_logger_st("debug_file", outputName + ".log");
    debugFileLogger->set_level(spdlog::level::debug);
  }
  debugFileLogger->debug("Debugging logger created");
  consoleLogger = spdlog::get("console");
  if (!consoleLogger) {
    consoleLogger = spdlog::stdout_color_st("console");
    consoleLogger->set_pattern("%v");
    consoleLogger->set_level(spdlog::level::warn);
  }
  debugFileLogger->debug("Console logger created");
  rawSpectrumLogger = spdlog::get("BSG_raw");
  if (!rawSpectrumLogger) {
    rawSpectrumLogger = spdlog::basic_logger_st("BSG_raw", outputName + ".raw");
    rawSpectrumLogger->set_pattern("%v");
    rawSpectrumLogger->set_level(spdlog::level::info);
  }
  debugFileLogger->debug("Raw spectrum logger created");
  resultsFileLogger = spdlog::get("BSG_results_file");
  if (!resultsFileLogger) {
    resultsFileLogger = spdlog::basic_logger_st("BSG_results_file", outputName + ".txt");
    resultsFileLogger->set_pattern("%v");
    resultsFileLogger->set_level(spdlog::level::info);
  }
  debugFileLogger->debug("Results file logger created");
}

void Generator::InitializeConstants() {
  debugFileLogger->debug("Entered initialize constants");

  Z = GetBSGOpt(int, Daughter.Z);
  A = GetBSGOpt(int, Daughter.A);

  R = GetBSGOpt(double, Daughter.Radius) * 1e-15 / NATURAL_LENGTH * std::sqrt(5. / 3.);
  if (R == 0.0) {
    debugFileLogger->debug("Radius not found. Using standard formula.");
    R = 1.2 * std::pow(A, 1. / 3.) * 1e-15 / NATURAL_LENGTH;
  }
  motherBeta2 = GetBSGOpt(double, Mother.Beta2);
  daughterBeta2 = GetBSGOpt(double, Daughter.Beta2);
  motherSpinParity = GetBSGOpt(int, Mother.SpinParity);
  daughterSpinParity = GetBSGOpt(int, Daughter.SpinParity);

  motherExcitationEn = GetBSGOpt(double, Mother.ExcitationEnergy);
  daughterExcitationEn = GetBSGOpt(double, Daughter.ExcitationEnergy);

  gA = GetBSGOpt(double, Constants.gA);
  gP = GetBSGOpt(double, Constants.gP);
  gM = GetBSGOpt(double, Constants.gM);

  std::string process = GetBSGOpt(std::string, Transition.Process);
  std::string type = GetBSGOpt(std::string, Transition.Type);

  if (boost::iequals(process, "B+")) {
    betaType = BETA_PLUS;
  } else {
    betaType = BETA_MINUS;
  }

  mixingRatio = 0.;
  if (boost::iequals(type, "Fermi")) {
    decayType = FERMI;
  } else if (boost::iequals(type, "Gamow-Teller")) {
    decayType = GAMOW_TELLER;
  } else {
    decayType = MIXED;
    mixingRatio = GetBSGOpt(double, Transition.MixingRatio);
  }

  if (A != GetBSGOpt(int, Mother.A)) {
    consoleLogger->error("Mother and daughter mass numbers are not the same.");
  }
  if (Z != GetBSGOpt(int, Mother.Z)+betaType) {
    consoleLogger->error("Mother and daughter cannot be obtained through {} process", process);
  }

  QValue = GetBSGOpt(double, Transition.QValue);

  atomicEnergyDeficit = GetBSGOpt(double, Transition.AtomicEnergyDeficit);

  if (betaType == BETA_MINUS) {
    W0 = (QValue - atomicEnergyDeficit + motherExcitationEn - daughterExcitationEn) / ELECTRON_MASS_KEV + 1.;
  } else {
    W0 = (QValue - atomicEnergyDeficit + motherExcitationEn - daughterExcitationEn) / ELECTRON_MASS_KEV - 1.;
  }
  W0 = W0 - (W0 * W0 - 1) / 2. / A / NUCLEON_MASS_KEV * ELECTRON_MASS_KEV;
  debugFileLogger->debug("Leaving InitializeConstants");
}

void Generator::InitializeShapeParameters() {
  debugFileLogger->debug("Entered InitializeShapeParameters");
  if (!OptExists(Spectrum.ModGaussFit)) {
    hoFit = CD::FitHODist(Z, R * std::sqrt(3. / 5.));
  } else {
    hoFit = GetBSGOpt(double, Spectrum.ModGaussFit);
  }
  debugFileLogger->debug("hoFit: {}", hoFit);

  baseShape = GetBSGOpt(std::string, Spectrum.Shape);

  vOld.resize(3);
  vNew.resize(3);

  if (baseShape == "Modified_Gaussian") {
    debugFileLogger->debug("Found Modified_Gaussian shape");
    vOld[0] = 3./2.;
    vOld[1] = -1./2.;
    vNew[0] = std::sqrt(5./2.)*4.*(1.+hoFit)*std::sqrt(2.+5.*hoFit)/std::sqrt(M_PI)*std::pow(2.+3.*hoFit, 3./2.);
    vNew[1] = -4./3./(3.*hoFit+2)/std::sqrt(M_PI)*std::pow(5.*(2.+5.*hoFit)/2./(2.+3.*hoFit), 3./2.);
    vNew[2] = (2.-7.*hoFit)/5./(3.*hoFit+2)/std::sqrt(M_PI)*std::pow(5.*(2.+5.*hoFit)/2./(2.+3.*hoFit), 5./3.);
  } else {
    if (OptExists(Spectrum.vold) && OptExists(Spectrum.vnew)) {
      debugFileLogger->debug("Found v and v'");
      vOld = GetBSGOpt(std::vector<double>, Spectrum.vold);
      vNew = GetBSGOpt(std::vector<double>, Spectrum.vnew);
    } else if (OptExists(vold) || OptExists(vnew)) {
      consoleLogger->error("ERROR: Both old and new potential expansions must be given.");
    }
  }
  debugFileLogger->debug("Leaving InitializeShapeParameters");
}

void Generator::LoadExchangeParameters() {
  debugFileLogger->debug("Entered LoadExchangeParameters");
  std::string exParamFile = GetBSGOpt(std::string, exchangedata);
  std::ifstream paramStream(exParamFile.c_str());
  std::string line;

  if (paramStream.is_open()) {
    while (getline(paramStream, line)) {
      double z;
      double a, b, c, d, e, f, g, h, i;

      std::istringstream iss(line);
      iss >> z >> a >> b >> c >> d >> e >> f >> g >> h >> i;

      if (z == Z - betaType) {
        exPars[0] = a;
        exPars[1] = b;
        exPars[2] = c;
        exPars[3] = d;
        exPars[4] = e;
        exPars[5] = f;
        exPars[6] = g;
        exPars[7] = h;
        exPars[8] = i;
      }
    }
  } else {
    consoleLogger->error("ERROR: Can't find Exchange parameters file at {}.", exParamFile);
  }
  debugFileLogger->debug("Leaving LoadExchangeParameters");
}

void Generator::InitializeL0Constants() {
  debugFileLogger->debug("Leaving InitializeL0Constants");
  double b[7][6];
  double bNeg[7][6];
  bNeg[0][0] = 0.115;
  bNeg[0][1] = -1.8123;
  bNeg[0][2] = 8.2498;
  bNeg[0][3] = -11.223;
  bNeg[0][4] = -14.854;
  bNeg[0][5] = 32.086;
  bNeg[1][0] = -0.00062;
  bNeg[1][1] = 0.007165;
  bNeg[1][2] = 0.01841;
  bNeg[1][3] = -0.53736;
  bNeg[1][4] = 1.2691;
  bNeg[1][5] = -1.5467;
  bNeg[2][0] = 0.02482;
  bNeg[2][1] = -0.5975;
  bNeg[2][2] = 4.84199;
  bNeg[2][3] = -15.3374;
  bNeg[2][4] = 23.9774;
  bNeg[2][5] = -12.6534;
  bNeg[3][0] = -0.14038;
  bNeg[3][1] = 3.64953;
  bNeg[3][2] = -38.8143;
  bNeg[3][3] = 172.1368;
  bNeg[3][4] = -346.708;
  bNeg[3][5] = 288.7873;
  bNeg[4][0] = 0.008152;
  bNeg[4][1] = -1.15664;
  bNeg[4][2] = 49.9663;
  bNeg[4][3] = -273.711;
  bNeg[4][4] = 657.6292;
  bNeg[4][5] = -603.7033;
  bNeg[5][0] = 1.2145;
  bNeg[5][1] = -23.9931;
  bNeg[5][2] = 149.9718;
  bNeg[5][3] = -471.2985;
  bNeg[5][4] = 662.1909;
  bNeg[5][5] = -305.6804;
  bNeg[6][0] = -1.5632;
  bNeg[6][1] = 33.4192;
  bNeg[6][2] = -255.1333;
  bNeg[6][3] = 938.5297;
  bNeg[6][4] = -1641.2845;
  bNeg[6][5] = 1095.358;

  double bPos[7][6];
  bPos[0][0] = 0.0701;
  bPos[0][1] = -2.572;
  bPos[0][2] = 27.5971;
  bPos[0][3] = -128.658;
  bPos[0][4] = 272.264;
  bPos[0][5] = -214.925;
  bPos[1][0] = -0.002308;
  bPos[1][1] = 0.066463;
  bPos[1][2] = -0.6407;
  bPos[1][3] = 2.63606;
  bPos[1][4] = -5.6317;
  bPos[1][5] = 4.0011;
  bPos[2][0] = 0.07936;
  bPos[2][1] = -2.09284;
  bPos[2][2] = 18.45462;
  bPos[2][3] = -80.9375;
  bPos[2][4] = 160.8384;
  bPos[2][5] = -124.8927;
  bPos[3][0] = -0.93832;
  bPos[3][1] = 22.02513;
  bPos[3][2] = -197.00221;
  bPos[3][3] = 807.1878;
  bPos[3][4] = -1566.6077;
  bPos[3][5] = 1156.3287;
  bPos[4][0] = 4.276181;
  bPos[4][1] = -96.82411;
  bPos[4][2] = 835.26505;
  bPos[4][3] = -3355.8441;
  bPos[4][4] = 6411.3255;
  bPos[4][5] = -4681.573;
  bPos[5][0] = -8.2135;
  bPos[5][1] = 179.0862;
  bPos[5][2] = -1492.1295;
  bPos[5][3] = 5872.5362;
  bPos[5][4] = -11038.7299;
  bPos[5][5] = 7963.4701;
  bPos[6][0] = 5.4583;
  bPos[6][1] = -115.8922;
  bPos[6][2] = 940.8305;
  bPos[6][3] = -3633.9181;
  bPos[6][4] = 6727.6296;
  bPos[6][5] = -4795.0481;

  for (int i = 0; i < 7; i++) {
    aPos[i] = 0;
    aNeg[i] = 0;
    for (int j = 0; j < 6; j++) {
      aNeg[i] += bNeg[i][j] * std::pow(ALPHA * Z, j + 1);
      aPos[i] += bPos[i][j] * std::pow(ALPHA * Z, j + 1);
    }
  }
  debugFileLogger->debug("Leaving InitializeL0Constants");
}

void Generator::InitializeNSMInfo() {
  nsm = new NS::NuclearStructureManager();

  if (OptExists(connect)) {
    int dKi, dKf;
    nsm->GetESPStates(spsi, spsf, dKi, dKf);
  }

  GetMatrixElements();
}

void Generator::GetMatrixElements() {
  debugFileLogger->info("Calculating matrix elements");
  double M101 = 1.0;
  if (!OptExists(Spectrum.Lambda)) {
    M101 = nsm->CalculateMatrixElement(false, 1, 0, 1);
    double M121 = nsm->CalculateMatrixElement(false, 1, 2, 1);
    ratioM121 = M121 / M101;
  } else {
    ratioM121 = GetBSGOpt(double, Spectrum.Lambda);
  }

  fc1 = gA * M101;

  bAc = dAc = 0;
  if (!OptExists(Spectrum.WeakMagnetism)) {
    debugFileLogger->info("Calculating Weak Magnetism");
    bAc = nsm->CalculateWeakMagnetism();
  } else {
    bAc = GetBSGOpt(double, Spectrum.WeakMagnetism);
  }
  if (!OptExists(Spectrum.InducedTensor)) {
    debugFileLogger->info("Calculating Induced Tensor");
    dAc = nsm->CalculateInducedTensor();
  } else {
    dAc = GetBSGOpt(double, Spectrum.InducedTensor);
  }
  debugFileLogger->info("Weak magnetism: {}", bAc);
  debugFileLogger->info("Induced tensor: {}", dAc);
  debugFileLogger->info("M121/M101: {}", ratioM121);

  fb = bAc * A * fc1;
  fd = dAc * A * fc1;
}

std::tuple<double, double> Generator::CalculateDecayRate(double W) {
  double result = 1;
  double neutrinoResult = 1;

  double Wv = W0 - W + 1;

  if (GetBSGOpt(bool, Spectrum.Phasespace)) {
    result *= SF::PhaseSpace(W, W0, motherSpinParity, daughterSpinParity);
    neutrinoResult *=
        SF::PhaseSpace(Wv, W0, motherSpinParity, daughterSpinParity);
  }

  if (GetBSGOpt(bool, Spectrum.Fermi)) {
    result *= SF::FermiFunction(W, Z, R, betaType);
    neutrinoResult *= SF::FermiFunction(Wv, Z, R, betaType);
  }
  if (GetBSGOpt(bool, Spectrum.C)) {
    if (OptExists(connect)) {
      result *= SF::CCorrection(W, W0, Z, A, R, betaType, hoFit, decayType, gA,
                                gP, fc1, fb, fd, ratioM121, GetBSGOpt(bool, Spectrum.Isovector), spsi, spsf);
      neutrinoResult *=
          SF::CCorrection(Wv, W0, Z, A, R, betaType, hoFit, decayType, gA, gP,
                          fc1, fb, fd, ratioM121, GetBSGOpt(bool, Spectrum.Isovector), spsi, spsf);
    } else {
      result *= SF::CCorrection(W, W0, Z, A, R, betaType, hoFit, decayType, gA,
                                gP, fc1, fb, fd, ratioM121, GetBSGOpt(bool, Spectrum.Isovector));
      neutrinoResult *=
          SF::CCorrection(Wv, W0, Z, A, R, betaType, hoFit, decayType, gA, gP,
                          fc1, fb, fd, ratioM121, GetBSGOpt(bool, Spectrum.Isovector));
    }
  }
  if (GetBSGOpt(bool, Spectrum.Relativistic)) {
    result *= SF::RelativisticCorrection(W, W0, Z, A, R, betaType, decayType);
    neutrinoResult *=
        SF::RelativisticCorrection(Wv, W0, Z, A, R, betaType, decayType);
  }
  if (GetBSGOpt(bool, Spectrum.ESDeformation)) {
    result *=
        SF::DeformationCorrection(W, W0, Z, R, daughterBeta2, betaType);
    neutrinoResult *=
        SF::DeformationCorrection(Wv, W0, Z, R, daughterBeta2, betaType);
  }
  if (GetBSGOpt(bool, Spectrum.ESFiniteSize)) {
    result *= SF::L0Correction(W, Z, R, betaType, aPos, aNeg);
    neutrinoResult *= SF::L0Correction(Wv, Z, R, betaType, aPos, aNeg);
  }
  if (GetBSGOpt(bool, Spectrum.ESShape)) {
    result *= SF::UCorrection(W, Z, R, betaType, baseShape, vOld, vNew);
    neutrinoResult *= SF::UCorrection(Wv, Z, R, betaType, baseShape, vOld, vNew);
  }
  if (GetBSGOpt(bool, Spectrum.CoulombRecoil)) {
    result *= SF::QCorrection(W, W0, Z, A, betaType, decayType, mixingRatio);
    neutrinoResult *=
        SF::QCorrection(Wv, W0, Z, A, betaType, decayType, mixingRatio);
  }
  if (GetBSGOpt(bool, Spectrum.Radiative)) {
    result *= SF::RadiativeCorrection(W, W0, Z, R, betaType, gA, gM);
    neutrinoResult *= SF::NeutrinoRadiativeCorrection(Wv);
  }
  if (GetBSGOpt(bool, Spectrum.Recoil)) {
    result *= SF::RecoilCorrection(W, W0, A, decayType, mixingRatio);
    neutrinoResult *= SF::RecoilCorrection(Wv, W0, A, decayType, mixingRatio);
  }
  if (GetBSGOpt(bool, Spectrum.Screening)) {
    result *= SF::AtomicScreeningCorrection(W, Z, betaType);
    neutrinoResult *= SF::AtomicScreeningCorrection(Wv, Z, betaType);
  }
  if (GetBSGOpt(bool, Spectrum.Exchange)) {
    if (betaType == BETA_MINUS) {
      result *= SF::AtomicExchangeCorrection(W, exPars);
      neutrinoResult *= SF::AtomicExchangeCorrection(Wv, exPars);
    }
  }
  if (GetBSGOpt(bool, Spectrum.AtomicMismatch)) {
    if (atomicEnergyDeficit == 0.) {
      result *= SF::AtomicMismatchCorrection(W, W0, Z, A, betaType);
      neutrinoResult *= SF::AtomicMismatchCorrection(Wv, W0, Z, A, betaType);
    }
  }
  result = std::max(0., result);
  neutrinoResult = std::max(0., neutrinoResult);
  rawSpectrumLogger->info("{:<10f}\t{:<10f}\t{:<10f}\t{:<10f}", W, (W-1.)*ELECTRON_MASS_KEV, result, neutrinoResult);

  return std::make_tuple(result, neutrinoResult);
}

std::vector<std::vector<double> > Generator::CalculateSpectrum() {
  debugFileLogger->info("Calculating spectrum");
  double beginEn = GetBSGOpt(double, Spectrum.Begin);
  double endEn = GetBSGOpt(double, Spectrum.End);
  double stepEn = GetBSGOpt(double, Spectrum.Step);

  double beginW = beginEn / ELECTRON_MASS_KEV + 1.;
  double endW = endEn / ELECTRON_MASS_KEV + 1.;
  if (endEn == 0.0) {
    endW = W0;
  }
  double stepW = stepEn / ELECTRON_MASS_KEV;

  double currentW = beginW;
  while (currentW <= endW) {
    auto result = CalculateDecayRate(currentW);
    std::vector<double> entry = {currentW, std::get<0>(result), std::get<1>(result)};
    spectrum.push_back(entry);
    currentW += stepW;
  }
  PrepareOutputFile();
  return spectrum;
}

double Generator::CalculateLogFtValue(double partialHalflife) {
  debugFileLogger->debug("Calculating Ft value with partial halflife {}", partialHalflife);
  double f = utilities::Simpson(spectrum);
  debugFileLogger->debug("f: {}", f);
  double ft = f*partialHalflife;
  double logFt = std::log10(ft);

  return logFt;
}

double Generator::CalculateMeanEnergy() {
  debugFileLogger->debug("Calculating mean energy");
  std::vector<std::vector<double> > weightedSpectrum;
  for (int i =0; i < spectrum.size(); i++) {
    std::vector<double> entry = {spectrum[i][0], spectrum[i][0]*spectrum[i][1]};
    weightedSpectrum.push_back(entry);
  }
  double weightedF = utilities::Simpson(weightedSpectrum);
  double f = utilities::Simpson(spectrum);
  debugFileLogger->debug("Weighted f: {} Clean f: {}", weightedF, f);
  return weightedF/f;
}

void Generator::PrepareOutputFile() {
  ShowBSGInfo();

  auto l = spdlog::get("BSG_results_file");
  l->info("Spectrum input overview\n{:=>30}", "");
  l->info("Using information from {}\n\n", GetBSGOpt(std::string, input));
  l->info("Transition from {}{} [{}/2] ({} keV) to {}{} [{}/2] ({} keV)", A, utilities::atoms[int(Z-1-betaType)], motherSpinParity, motherExcitationEn, A, utilities::atoms[int(Z-1)], daughterSpinParity, daughterExcitationEn);
  l->info("Q Value: {} keV\tEffective endpoint energy: {}", QValue, (W0-1.)*ELECTRON_MASS_KEV);
  l->info("Process: {}\tType: {}", GetBSGOpt(std::string, Transition.Process), GetBSGOpt(std::string, Transition.Type));
  if (mixingRatio != 0) l->info("Mixing ratio: {}", mixingRatio);
  if (OptExists(Transition.PartialHalflife)) {
    l->info("Partial halflife: {} s", GetBSGOpt(double, Transition.PartialHalflife));
    l->info("Calculated log ft value: {}", CalculateLogFtValue(GetBSGOpt(double, Transition.PartialHalflife)));
  } else {
    l->info("Partial halflife: not given");
    l->info("Calculated log f value: {}", CalculateLogFtValue(1.0));
  }
  l->info("Mean energy: {} keV", (CalculateMeanEnergy()-1.)*ELECTRON_MASS_KEV);
  l->info("\nMatrix Element Summary\n{:->30}", "");
  if (OptExists(Spectrum.WeakMagnetism)) l->info("{:25}: {} ({})", "b/Ac (weak magnetism)", bAc, "given");
  else l->info("{:35}: {}", "b/Ac (weak magnetism)", bAc);
  if (OptExists(Spectrum.Inducedtensor)) l->info("{:25}: {} ({})", "d/Ac (induced tensor)", dAc, "given");
  else l->info("{:35}: {}", "d/Ac (induced tensor)", dAc);
  if (OptExists(Spectrum.Lambda)) l->info("{:25}: {} ({})", "AM121/AM101", ratioM121, "given");
  else l->info("{:35}: {}", "AM121/AM101", ratioM121);

  l->info("Full breakfown written in {}.nme", GetBSGOpt(std::string, output));

  l->info("\nSpectral corrections\n{:->30}", "");
  l->info("{:25}: {}", "Phase space", GetBSGOpt(bool, Spectrum.Phasespace));
  l->info("{:25}: {}", "Fermi function", GetBSGOpt(bool, Spectrum.Fermi));
  l->info("{:25}: {}", "L0 correction", GetBSGOpt(bool, Spectrum.ESFiniteSize));
  l->info("{:25}: {}", "C correction", GetBSGOpt(bool, Spectrum.C));
  l->info("{:25}: {}", "Isovector correction", GetBSGOpt(bool, Spectrum.Isovector));
  l->info("    Connected: {}", GetBSGOpt(bool, Spectrum.Connect));
  l->info("{:25}: {}", "Relativistic terms", GetBSGOpt(bool, Spectrum.Relativistic));
  l->info("{:25}: {}", "Deformation", GetBSGOpt(bool, Spectrum.ESDeformation));
  l->info("{:25}: {}", "U correction", GetBSGOpt(bool, Spectrum.ESShape));
  l->info("    Shape: {}", GetBSGOpt(std::string, Spectrum.Shape));
  if (OptExists(Spectrum.vold) && OptExists(Spectrum.vnew)) {
    l->info("    v : {}, {}, {}", vOld[0], vOld[1], vOld[2]);
    l->info("    v': {}, {}, {}", vNew[0], vNew[1], vNew[2]);
  } else {
    l->info("    v : not given");
    l->info("    v': not given");
  }
  l->info("{:25}: {}", "Q correction", GetBSGOpt(bool, Spectrum.CoulombRecoil));
  l->info("{:25}: {}", "Radiative correction", GetBSGOpt(bool, Spectrum.Radiative));
  l->info("{:25}: {}", "Nuclear recoil", GetBSGOpt(bool, Spectrum.Recoil));
  l->info("{:25}: {}", "Atomic screening", GetBSGOpt(bool, Spectrum.Screening));
  l->info("{:25}: {}", "Atomic exchange", GetBSGOpt(bool, Spectrum.Exchange));
  l->info("{:25}: {}", "Atomic mismatch", GetBSGOpt(bool, Spectrum.AtomicMismatch));
  l->info("{:25}: {}", "Export neutrino", GetBSGOpt(bool, Spectrum.Neutrino));

  l->info("\n\nSpectrum calculated from {} keV to {} keV with step size {} keV\n", GetBSGOpt(double, Spectrum.Begin), GetBSGOpt(double, Spectrum.End) > 0 ? GetBSGOpt(double, Spectrum.End) : (W0-1.)*ELECTRON_MASS_KEV, GetBSGOpt(double, Spectrum.Step));

  if (GetBSGOpt(bool, Spectrum.Neutrino))  l->info("{:10}\t{:10}\t{:10}\t{:10}", "W [m_ec2]", "E [keV]", "dN_e/dW", "dN_v/dW");
  else l->info("{:10}\t{:10}\t{:10}", "W [m_ec2]", "E [keV]", "dN_e/dW");

  for (int i = 0; i < spectrum.size(); i++) {
    if (GetBSGOpt(bool, Spectrum.Neutrino)) {
      l->info("{:<10f}\t{:<10f}\t{:<10f}\t{:<10f}", spectrum[i][0], (spectrum[i][0]-1.)*ELECTRON_MASS_KEV, spectrum[i][1], spectrum[i][2]);
    } else {
      l->info("{:<10f}\t{:<10f}\t{:<10f}", spectrum[i][0], (spectrum[i][0]-1.)*ELECTRON_MASS_KEV, spectrum[i][1]);
    }
  }
}
