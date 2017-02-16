#include "OptionContainer.h"
#include "Constants.h"
#include "boost/algorithm/string.hpp"

#define GetOpt(a, b) OptionContainer::GetInstance().GetOption< a >( #b )

Generator::Generator() {
  Z = GetOpt(int, NuclearProperties.DaughterZ);
  A = GetOpt(int, NuclearProperties.DaughterA);
  R = GetOpt(double, NuclearProperties.DaughterRadius)*1e-15/NATLENGTH;
  motherBeta2 = GetOpt(double, NuclearProperties.MotherBeta2);
  motherBeta4 = GetOpt(double, NuclearProperties.MotherBeta4);
  daughterBeta2 = GetOpt(double, NuclearProperties.DaughterBeta2);
  daughterBeta4 = GetOpt(double, NuclearProperties.DaughterBeta4);
  motherSpinParity = GetOpt(std::string, NuclearProperties.MotherSpinParity);
  daughterSpinParity = GetOpt(std::string, NuclearProperties.DaughterSpinParity);
  
  gA = GetOpt(double, Constants.gA);
  gP = GetOpt(double, Constants.gP);
  gM = GetOpt(double, Constants.gM);

  std::string process = getOpt(std::string, Transition.Process);
  std::string type = GetOpt(std::string, Transition.Type);

  if (boost::iequals(process, "B+")) {
    fBetaType = BETA_PLUS;
  }
  else {
    fBetaType = BETA_MINUS;
  }

  mixingRatio = 0.;
  if (boost::iequals(type, "fermi")) {
    fDecayType = FERMI;
  }
  else if (boost::iequals(type, "gamow-teller")) {
    fDecayType = GAMOW_TELLER;
  }
  else 
    fDecayType = MIXED;
    mixingRatio = GetOpt(double, Transition.MixingRatio);
  }

  double QValue = GetOpt(double, Transition.QValue);
  
  if (fBetaType == BETA_MINUS) {
    W0 = QValue/electronMasskeV + 1.;
  } 
  else {
    W0 = QValue/electronMasskeV - 1.;
  }
  W0 = W0 - (W0*W0-1)/2./A/1837.4;

  LoadExchangeParameters();

  InitializeL0Constants();

  CalculateMatrixElements();
}

void Generator::LoadExchangeParameters() {
  exParamFile = GetOpt(std::string, ExchangeData);
  ifstream paramStream(exParamFile.c_str());
  string line;

  if (paramStream.is_open()) {
    while (getline(paramStream, line)) {
      double Z;
      double a, b, c, d, e, f, g, h, i;

      istringstream iss(line);
      iss >> Z >> a >> b >> c >> d >> e >> f >> g >> h >> i;

      if (Z == Zd - fBetaType) {
        exPars[0] = a;
        exPars[1] = b;
        exPars[2] = c;
        exPars[3] = d;
        exPars[4] = e;
        exPars[5] = f;
        exPars[6] = g;
        exPars[7] = h;
        exPars[8] = i;

        // cout << "Loaded Exchange Params for Z=" << Z << ": " << a << " " << b
        //     << " " << c << " and so on" << endl;
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
  utilities::NuclearState nsf = nilsson::CalculateDeformedState(Z, A, daughterBeta2, daughterBeta4);
  utilities::NuclearState nsi =
      nilsson::CalculateDeformedState(Z + fBetaType, A, motherBeta2, motherBeta4);

  return -std::sqrt(2. / 3.) * (protonMasskeV + neutronMasskeV) / 2. /
             electronMasskeV * R / gA *
             utilities::CalculateDeformedSPMatrixElement(
                 nsi.states, nsf.states, true, 1, 1, 1, nsi.O, nsf.O, nsi.K,
                 nsf.K, R) /
             utilities::CalculateDeformedSPMatrixElement(
                 nsi.states, nsf.states, false, 1, 0, 1, nsi.O, nsf.O, nsi.K,
                 nsf.K, R) +
         gM / gA;
}

double Generator::CalculateMatrixElements() {
  
}

vector<double> SpectralFunctions::GetSpectrumShape(double ) {
  double W = Energy / electronMasskeV + 1;

  double result = 1;
  double neutrinoResult = 1;

  double Wv = W0 - W + 1;

  if (fPhaseSpace) {
    result *= std::std::sqrt(W * W - 1) * W * sqr(W0 - W);
    neutrinoResult *= std::std::sqrt(Wv * Wv - 1) * Wv * sqr(W0 - Wv);
  }

  if (fFermiFunction) {
    result *= FermiFunction(W);
    neutrinoResult *= FermiFunction(Wv);
  }

  if (fCCorrection) {
    result *= CCorrection(W);
    neutrinoResult *= CCorrection(Wv);
  }

  /*if (fQCDInducedCorrection) {
    result *= QCDInducedCorrection(W);
    neutrinoResult *= QCDInducedCorrection(Wv);
  }*/

  if (fRelativisticCorrection) {
    result *= RelativisticCorrection(W);
    neutrinoResult *= RelativisticCorrection(Wv);
  }

  if (fDeformationCorrection) {
    result *= DeformationCorrection(W);
    neutrinoResult *= DeformationCorrection(Wv);
  }

  if (fL0Correction) {
    result *= L0Correction(W, R);
    neutrinoResult *= L0Correction(Wv, R);
  }

  if (fUCorrection) {
    result *= UCorrection(W);
    neutrinoResult *= UCorrection(Wv);
  }

  if (fQCorrection) {
    result *= QCorrection(W);
    neutrinoResult *= QCorrection(Wv);
  }

  if (fRadiativeCorrection) {
    result *= RadiativeCorrection(W);
    neutrinoResult *= NeutrinoRadiativeCorrection(Wv);
  }

  if (fRecoilCorrection) {
    result *= RecoilCorrection(W);
    neutrinoResult *= RecoilCorrection(Wv);
  }

  if (fAtomicScreeningCorrection) {
    result *= AtomicScreeningCorrection(W);
    neutrinoResult *= AtomicScreeningCorrection(Wv);
  }

  if (fAtomicExchAgeCorrection) {
    if (fBetaType == BETA_MINUS) {
      result *= AtomicExchAgeCorrection(W);
      neutrinoResult *= AtomicExchAgeCorrection(Wv);
    } else {
      // cout << "The exchAge correction does not exist for beta+ decay." <<
      // endl;
    }
  }

  if (fAtomicMismatchCorrection) {
    result *= AtomicMismatchCorrection(W);
    neutrinoResult *= AtomicMismatchCorrection(Wv);
  }

  vector<double> fullResult;
  fullResult.push_back(result);
  fullResult.push_back(neutrinoResult);
  return fullResult;
}

void SpectralFunctions::WriteSpectrumToStream(ostream *ofile, double Step) {
  double CurrEn = Step;
  ofile->setf(std::ios::scientific);
  do {
    vector<double> spectShape = GetSpectrumShape(CurrEn);
    if (fCalcNeutrinoSpectrum)
      *ofile << CurrEn << "\t" << spectShape[0] << "\t" << spectShape[1]
             << endl;
    else
      *ofile << CurrEn << "\t" << spectShape[0] << endl;
    CurrEn += Step;
  } while (CurrEn < EndPointEnergy);
  return;
}
