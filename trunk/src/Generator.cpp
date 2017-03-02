#include "Generator.h"

#include "OptionContainer.h"
#include "ChargeDistributions.h"
#include "Constants.h"
#include "Utilities.h"
#include "MatrixElements.h"
#include "SpectralFunctions.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "boost/algorithm/string.hpp"

namespace SF = SpectralFunctions;
namespace ME = MatrixElements;
namespace CD = ChargeDistributions;

using std::cout;
using std::endl;

Generator::Generator() {
  Z = GetOpt(int, NuclearProperties.DaughterZ);
  A = GetOpt(int, NuclearProperties.DaughterA);
  R = GetOpt(double, NuclearProperties.DaughterRadius) * 1e-15 / NATLENGTH * std::sqrt(5./3.); 
  if (R == 0.0) {
    R = 1.2 * std::pow(A, 1./3.) * 1e-15 / NATLENGTH;
  }
  motherBeta2 = GetOpt(double, NuclearProperties.MotherBeta2);
  motherBeta4 = GetOpt(double, NuclearProperties.MotherBeta4);
  daughterBeta2 = GetOpt(double, NuclearProperties.DaughterBeta2);
  daughterBeta4 = GetOpt(double, NuclearProperties.DaughterBeta4);
  motherSpinParity = GetOpt(int, NuclearProperties.MotherSpinParity);
  daughterSpinParity = GetOpt(int, NuclearProperties.DaughterSpinParity);

  gA = GetOpt(double, Constants.gA);
  gP = GetOpt(double, Constants.gP);
  gM = GetOpt(double, Constants.gM);

  std::string process = GetOpt(std::string, Transition.Process);
  std::string type = GetOpt(std::string, Transition.Type);

  if (boost::iequals(process, "B+")) {
    fBetaType = BETA_PLUS;
  } else {
    fBetaType = BETA_MINUS;
  }

  mixingRatio = 0.;
  if (boost::iequals(type, "Fermi")) {
    fDecayType = FERMI;
  } else if (boost::iequals(type, "Gamow-Teller")) {
    fDecayType = GAMOW_TELLER;
  } else {
    fDecayType = MIXED;
    mixingRatio = GetOpt(double, Transition.MixingRatio);
  }

  double QValue = GetOpt(double, Transition.QValue);

  if (fBetaType == BETA_MINUS) {
    W0 = QValue / electronMasskeV + 1.;
  } else {
    W0 = QValue / electronMasskeV - 1.;
  }
  W0 = W0 - (W0 * W0 - 1) / 2. / A / nucleonMasskeV * electronMasskeV;

  LoadExchangeParameters();

  InitializeL0Constants();

  CalculateMatrixElements();

  hoFit = CD::FitHODist(Z, R * std::sqrt(3. / 5.));
}

void Generator::LoadExchangeParameters() {
  std::string exParamFile = GetOpt(std::string, ExchangeData);
  std::ifstream paramStream(exParamFile.c_str());
  std::string line;

  if (paramStream.is_open()) {
    while (getline(paramStream, line)) {
      double Z;
      double a, b, c, d, e, f, g, h, i;

      std::istringstream iss(line);
      iss >> Z >> a >> b >> c >> d >> e >> f >> g >> h >> i;

      if (Z == Z - fBetaType) {
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
  }
}

void Generator::InitializeL0Constants() {
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
      aNeg[i] += bNeg[i][j] * std::pow(alpha * Z, j + 1);
      aPos[i] += bPos[i][j] * std::pow(alpha * Z, j + 1);
    }
  }
}

double Generator::CalculateWeakMagnetism() {
  double result = 0.0;
  std::string pot = GetOpt(std::string, Computational.Potential);
  if (boost::iequals(pot, "SHO")) {
    cout << "Weak magnetism SHO" << endl;
    int ni, li, si, nf, lf, sf;
    GetSPOrbitalNumbers(ni, li, si, nf, lf, sf);
    cout << li << " " << si << " " << lf << " " << sf << endl;
    if (li == lf) {
      if (si == sf && si > 0) {
        return 1. / gA * (li + 1 + gM);
      } else if (si == sf) {
        return -1. / gA * (li - gM);
      } else {
        return 1. / 2. / gA;
      }
    }
  } else {
    double nu = CD::CalcNu(R * std::sqrt(3./5.), Z);
    result =
        -std::sqrt(2. / 3.) * nucleonMasskeV / electronMasskeV * R / gA *
            ME::CalculateDeformedSPMatrixElement(
                spsi, spsf, true, 1, 1, 1, spsi.O, spsf.O, spsi.K, spsf.K, R, nu) /
            ME::CalculateDeformedSPMatrixElement(
                spsi, spsf, false, 1, 0, 1, spsi.O, spsf.O, spsi.K, spsf.K, R, nu) +
        gM / gA;
  }
  return result;
}

double Generator::CalculateInducedTensor() {
  double result = 0.0;
  std::string pot = GetOpt(std::string, Computational.Potential);
  double nu = CD::CalcNu(R * std::sqrt(3./5.), Z);
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
  } else {
    result =
        2. / std::sqrt(3.) * nucleonMasskeV / electronMasskeV * R *
        ME::CalculateDeformedSPMatrixElement(spsi, spsf, false, 1, 1, 0, spsi.O,
                                             spsf.O, spsi.K, spsf.K, R, nu) /
        ME::CalculateDeformedSPMatrixElement(spsi, spsf, false, 1, 0, 1, spsi.O,
                                             spsf.O, spsi.K, spsf.K, R, nu);
  }
  return result;
}

double Generator::CalculateRatioM121() {
  double M121, M101;
  std::string pot = GetOpt(std::string, Computational.Potential);
  cout << "Potential: " << pot << endl;
  double nu = CD::CalcNu(R * std::sqrt(3. / 5.), Z);
  if (boost::iequals(pot, "SHO")) {
    int ni, li, si, nf, lf, sf;
    GetSPOrbitalNumbers(ni, li, si, nf, lf, sf);
    M121 = ME::GetSingleParticleMatrixElement(
        false, std::abs(motherSpinParity) / 2., 1, 2, 1, ni, nf, li, lf, si, sf, R, nu);
    M101 = ME::GetSingleParticleMatrixElement(
        false, std::abs(motherSpinParity) / 2., 1, 0, 1, ni, nf, li, lf, si, sf, R, nu);

  } else {
    M121 = ME::CalculateDeformedSPMatrixElement(
        spsi, spsf, false, 1, 2, 1, spsi.O, spsf.O, spsi.K, spsf.K, R, nu);
    M101 = ME::CalculateDeformedSPMatrixElement(
        spsi, spsf, false, 1, 0, 1, spsi.O, spsf.O, spsi.K, spsf.K, R, nu);
  }
  cout << "M121: " << M121 << " M101: " << M101 << endl;
  fc1 = gA * M101;

  return M121 / M101;
}

void Generator::GetSPOrbitalNumbers(int& ni, int& li, int& si, int& nf, int& lf,
                                    int& sf) {
  std::vector<int> occNumbersInit, occNumbersFinal;
  if (fDecayType == BETA_MINUS) {
    occNumbersInit = utilities::GetOccupationNumbers(A - (Z - fBetaType));
    occNumbersFinal = utilities::GetOccupationNumbers(Z);
  } else {
    occNumbersInit = utilities::GetOccupationNumbers(Z - fBetaType);
    occNumbersFinal = utilities::GetOccupationNumbers(A - Z);
  }
  ni = occNumbersInit[occNumbersInit.size() - 1 - 3];
  li = occNumbersInit[occNumbersInit.size() - 1 - 2];
  si = occNumbersInit[occNumbersInit.size() - 1 - 1];
  nf = occNumbersFinal[occNumbersFinal.size() - 1 - 3];
  lf = occNumbersFinal[occNumbersFinal.size() - 1 - 2];
  sf = occNumbersFinal[occNumbersFinal.size() - 1 - 1];
}

void Generator::GetDeformedStates(nilsson::SingleParticleState& spsf,
                                  nilsson::SingleParticleState& spsi,
                                  std::string pot) {
  double V0p = GetOpt(double, Computational.V0proton);
  double V0n = GetOpt(double, Computational.V0neutron);
  double A0 = GetOpt(double, Computational.SurfaceThickness);
  double VSp = GetOpt(double, Computational.V0Sproton);
  double VSn = GetOpt(double, Computational.V0Sneutron);

  double r = R * NATLENGTH * 1e15;

  if (boost::iequals(pot, "WS")) {
    if (fDecayType == BETA_MINUS) {
      spsf =
          nilsson::CalculateDeformedState(Z, 0, A, r, 0.0, 0.0, V0p, A0, VSp);
      spsi = nilsson::CalculateDeformedState(0, A - (Z - fDecayType), A, r, 0.0,
                                             0.0, V0n, A0, VSn);
    } else {
      spsf = nilsson::CalculateDeformedState(0, A - Z, A, r, 0.0, 0.0, V0n, A0,
                                             VSn);
      spsi = nilsson::CalculateDeformedState(Z - fDecayType, 0, A, r, 0.0, 0.0,
                                             V0p, A0, VSp);
    }
  } else if (boost::iequals(pot, "DWS")) {
    if (fDecayType == BETA_MINUS) {
      spsf = nilsson::CalculateDeformedState(Z, 0, A, r, daughterBeta2,
                                             daughterBeta4, V0p, A0, VSp);
      spsi = nilsson::CalculateDeformedState(0, A - (Z - fDecayType), A, r,
                                             motherBeta2, motherBeta4, V0n, A0,
                                             VSn);
    } else {
      spsf = nilsson::CalculateDeformedState(0, A - Z, A, r, daughterBeta2,
                                             daughterBeta4, V0n, A0, VSn);
      spsi = nilsson::CalculateDeformedState(
          Z - fDecayType, 0, A, r, motherBeta2, motherBeta4, V0p, A0, VSp);
    }
  }
}

void Generator::CalculateMatrixElements() {
  cout << "Calculating matrix elements" << endl;
  GetDeformedStates(spsf, spsi, GetOpt(std::string, Computational.Potential));
  ratioM121 = CalculateRatioM121();

  cout << "Ratio: " << ratioM121 << endl;

  double bAc = 0, dAc = 0;
  if (!OptExists(weakmagnetism)) {
    cout << "Calculating Weak Magnetism" << endl;
    bAc = CalculateWeakMagnetism();
  }
  if (!OptExists(inducedtensor)) {
    dAc = CalculateInducedTensor();
  }
  cout << "Weak magnetism: " << bAc << endl;
  cout << "Induced tensor: " << dAc << endl;
  cout << "M121/M101: " << ratioM121 << endl;

  fb = bAc * A * fc1;
  fd = dAc * A * fc1;
}

double* Generator::CalculateDecayRate(double W) {
  double result = 1;
  double neutrinoResult = 1;

  double Wv = W0 - W + 1;

  if (!OptExists(phase)) {
    result *= SF::PhaseSpace(W, W0, motherSpinParity, daughterSpinParity);
    neutrinoResult *=
        SF::PhaseSpace(Wv, W0, motherSpinParity, daughterSpinParity);
  }

  cout << "result: " << result << endl;

  if (!OptExists(fermi)) {
    result *= SF::FermiFunction(W, Z, R, fBetaType);
    neutrinoResult *= SF::FermiFunction(Wv, Z, R, fBetaType);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(C)) {
    result *= SF::CCorrection(W, W0, Z, A, R, fBetaType, hoFit, fDecayType, gA,
                              gP, fc1, fb, fd, ratioM121);
    neutrinoResult *=
        SF::CCorrection(Wv, W0, Z, A, R, fBetaType, hoFit, fDecayType, gA, gP,
                        fc1, fb, fd, ratioM121);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(relativistic)) {
    result *= SF::RelativisticCorrection(W, W0, Z, A, R, fBetaType, fDecayType);
    neutrinoResult *=
        SF::RelativisticCorrection(Wv, W0, Z, A, R, fBetaType, fDecayType);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(deformation)) {
    result *=
        SF::DeformationCorrection(W, W0, Z, R, daughterBeta2, daughterBeta4);
    neutrinoResult *=
        SF::DeformationCorrection(Wv, W0, Z, R, daughterBeta2, daughterBeta4);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(l0)) {
    result *= SF::L0Correction(W, Z, R, fBetaType, aPos, aNeg);
    neutrinoResult *= SF::L0Correction(Wv, Z, R, fBetaType, aPos, aNeg);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(U)) {
    result *= SF::UCorrection(W, Z, fBetaType);
    neutrinoResult *= SF::UCorrection(Wv, Z, fBetaType);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(Q)) {
    result *= SF::QCorrection(W, W0, Z, A, fBetaType, fDecayType, mixingRatio);
    neutrinoResult *=
        SF::QCorrection(Wv, W0, Z, A, fBetaType, fDecayType, mixingRatio);
  }

   cout << "result: " << result << endl;

  if (!OptExists(radiative)) {
    result *= SF::RadiativeCorrection(W, W0, Z, R, fBetaType, gA, gM);
    neutrinoResult *= SF::NeutrinoRadiativeCorrection(Wv);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(recoil)) {
    result *= SF::RecoilCorrection(W, W0, A, fDecayType, mixingRatio);
    neutrinoResult *= SF::RecoilCorrection(Wv, W0, A, fDecayType, mixingRatio);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(screening)) {
    result *= SF::AtomicScreeningCorrection(W, Z, fBetaType);
    neutrinoResult *= SF::AtomicScreeningCorrection(Wv, Z, fBetaType);
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(exchange)) {
    if (fBetaType == BETA_MINUS) {
      result *= SF::AtomicExchangeCorrection(W, exPars);
      neutrinoResult *= SF::AtomicExchangeCorrection(Wv, exPars);
    }
  }
  cout << "result: " << result << endl;

 
  if (!OptExists(mismatch)) {
    result *= SF::AtomicMismatchCorrection(W, W0, Z, A, fBetaType);
    neutrinoResult *= SF::AtomicMismatchCorrection(Wv, W0, Z, A, fBetaType);
  }
  cout << "result: " << result << endl;

 
  double fullResult[2] = {result, neutrinoResult};
  return fullResult;
}

std::vector<std::vector<double> > Generator::CalculateSpectrum() {
  double beginEn = GetOpt(double, begin);
  double endEn = GetOpt(double, end);
  double stepEn = GetOpt(double, step);

  double beginW = beginEn / electronMasskeV + 1.;
  double endW = endEn / electronMasskeV + 1.;
  if (endEn == 0.0) {
    endW = W0;
  }
  double stepW = stepEn / electronMasskeV;

  double currentW = beginW;
  while (currentW <= endW) {
    double* result = CalculateDecayRate(currentW);
    std::vector<double> entry = {currentW, result[0], result[1]};
    spectrum.push_back(entry);
    currentW += stepW;
  }
  return spectrum;
}

void Generator::WriteSpectrumToFile() {
  std::string outputName = GetOpt(std::string, output);
  std::ofstream outputStream(outputName.c_str());

  outputStream.setf(std::ios::scientific);
  for (int i = 0; i < spectrum.size(); i++) {
    outputStream << spectrum[i][0] << "\t" << spectrum[i][1];
    if (!OptExists(neutrino)) {
      outputStream << "\t" << spectrum[i][2];
    }
    outputStream << endl;
  }
}
