Spectrum::Spectrum() {
  An = 60;
  Zd = 28;
  EndPointEnergy = 320;  // in keV
  MixingRatio = 0.;
  fBetaType = BETA_MINUS;
  fDecayType = GAMOW_TELLER;
  fPhaseSpace = true;
  fFermiFunction = true;
  fCCorrection = true;
  fCICorrection = true;
  // fQCDInducedCorrection = false;
  fRelativisticCorrection = true;
  fDeformationCorrection = true;
  fL0Correction = true;
  fUCorrection = true;
  fQCorrection = true;
  fRadiativeCorrection = true;
  fRecoilCorrection = true;
  fAtomicScreeningCorrection = true;
  fAtomicExchangeCorrection = true;
  fAtomicMismatchCorrection = true;

  Afit = -1;

  fCalcNeutrinoSpectrum = false;

  exParamFile = "data/ExchangeData.dat";

  fb = 0;
  fc1 = 1.;
  fd = 0;

  gA = 1.1;
  gP = 0.;
  gM = 4.706;

  beta2 = 0.;
}

void Spectrum::ReadFromFile(char *FileName) {
  ifstream ifile(FileName, ifstream::in);
  if (!ifile.is_open()) {
    cerr << "Error opening " << FileName << endl;
    return;
  }
  float he;  /// temporary variable for sscanf
  char fline[100];
  ifile.getline(fline, 100);  // skip description
  ifile.getline(fline, 100);
  sscanf(fline, "%f", &he);  // must read into float type
  An = he;
  ifile.getline(fline, 100);
  ifile.getline(fline, 100);
  sscanf(fline, "%f", &he);
  Zd = he;
  ifile.getline(fline, 100);
  ifile.getline(fline, 100);
  sscanf(fline, "%f", &he);
  double QValue = he;
  // decay type
  int temp;
  ifile.getline(fline, 100);
  ifile.getline(fline, 100);
  sscanf(fline, "%i", &temp);
  if (temp == 0) {
    fBetaType = BETA_MINUS;
  } else {
    fBetaType = BETA_PLUS;
  }
  ifile.getline(fline, 100);
  ifile.getline(fline, 100);
  sscanf(fline, "%i", &temp);
  if (temp == 0) {
    fDecayType = FERMI;
  } else if (temp == 1) {
    fDecayType = GAMOW_TELLER;
  } else {
    fDecayType = MIXED;
    ifile.getline(fline, 100);
    ifile.getline(fline, 100);
    float t;
    sscanf(fline, "%f", &t);
    MixingRatio = t;
  }
  ifile.getline(fline, 100);
  ifile.getline(fline, 100);
  sscanf(fline, "%f", &he);

  if (he != QValue) {
    R = sqrt(5. / 3.) * he * 1e-15 * electronMasskeV * 1000 / hbar /
        speed_of_light;
  } else {
    R = 1.2 * pow(An, 1. / 3.) * 1e-15 * electronMasskeV * 1000 / hbar /
        speed_of_light;
  }

  he = 0.;

  ifile.getline(fline, 100);
  ifile.getline(fline, 100);
  sscanf(fline, "%f", &he);

  beta2 = he;

  if (fBetaType == BETA_MINUS) {
    W0 = QValue / electronMasskeV + 1.;
  } else {
    W0 = QValue / electronMasskeV - 1.;
  }
  W0 = W0 - (W0 * W0 - 1) / 2. / An / 1837.4;
  EndPointEnergy = (W0 - 1.) * electronMasskeV;

  if (bAc != 0) {
    fc1 = 1.;
    fb = An * bAc;
  } else {
    fc1 = 1.;
    bAc = utilities::CalculateWeakMagnetism(gA, gM, R, Zd, An, beta2,
                                            fBetaType);
    cout << "Received weak magnetism " << bAc << endl;
    fb = An * bAc;
  }

  ifile.close();

  LoadExchangeParameters();
}

void Spectrum::LoadExchangeParameters() {
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

void Spectrum::Initialize() {
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
      aNeg[i] += bNeg[i][j] * pow(alpha * Zd, j + 1);
      aPos[i] += bPos[i][j] * pow(alpha * Zd, j + 1);
    }
  }

  W0 = EndPointEnergy / electronMasskeV + 1;
  gamma = sqrt(1 - sqr(alpha * Zd));

  /*double r_0 = (1.15 + 1.80 * pow(An, -2. / 3.) - 1.20 * pow(An, -4. / 3.)) *
               1e-15 * electronMasskeV * 1000 / hbar / speed_of_light;
  R = r_0 * pow(An, 1. / 3.);*/
}

void Spectrum::PrintStatus(std::ostream *ofile) {
  *ofile << "SpectrumGenerator Status Information: " << endl;
  *ofile << "\t\t A=" << An << endl;
  *ofile << "\t\t Z=" << Zd << endl;
  *ofile << "\t\t E0=" << EndPointEnergy << " keV" << endl;
  //*ofile << "\t\t Ntilde (screening)=" << Ntilde << endl;
  if (fBetaType == BETA_PLUS) {
    *ofile << "\t\t positive beta decay" << endl;
  } else {
    *ofile << "\t\t negative beta decay" << endl;
  }
  *ofile << "\t\t Nuclear radius R = " << R << " (natural units)" << endl;
  *ofile << "\t\t Deformation: " << beta2 << endl;
  if (fDecayType == FERMI) {
    *ofile << "\t\t Pure Fermi transition" << endl;
  } else if (fDecayType == GAMOW_TELLER) {
    *ofile << "\t\t Pure Gamow-Teller transition" << endl;
  } else {
    *ofile << "\t\t Mixed Fermi/Gamow-Teller transition. Mixing ratio: "
           << MixingRatio << endl;
  }

  *ofile << "\t Options:" << endl;
  *ofile << "\t\t Phase Space factor: " << (fPhaseSpace ? "true" : "false")
         << endl;
  *ofile << "\t\t Fermi Function: " << (fFermiFunction ? "true" : "false")
         << endl;
  *ofile << "\t\t C Correction: " << (fCCorrection ? "true" : "false") << endl;
  *ofile << "\t\t CI Correction: " << (fCICorrection ? "true" : "false")
         << endl;
  /**ofile << "\t\t QCD Induced Correction: "
         << (fQCDInducedCorrection ? "true" : "false") << endl;*/
  *ofile << "\t\t Relativistic Correction: "
         << (fRelativisticCorrection ? "true" : "false") << endl;
  *ofile << "\t\t Deformation Correction: "
         << (fDeformationCorrection ? "true" : "false") << endl;
  *ofile << "\t\t L0 Correction: " << (fL0Correction ? "true" : "false")
         << endl;
  *ofile << "\t\t U Correction: " << (fUCorrection ? "true" : "false") << endl;
  *ofile << "\t\t Q Correction: " << (fQCorrection ? "true" : "false") << endl;
  *ofile << "\t\t Radiative Correction: "
         << (fRadiativeCorrection ? "true" : "false") << endl;
  *ofile << "\t\t Recoil Correction: " << (fRecoilCorrection ? "true" : "false")
         << endl;
  *ofile << "\t\t Atomic Screening: "
         << (fAtomicScreeningCorrection ? "true" : "false") << endl;
  *ofile << "\t\t Atomic Exchange Correction: "
         << (fAtomicExchangeCorrection ? "true" : "false") << endl;
  *ofile << "\t\t Atomic Mismatch Correction: "
         << (fAtomicMismatchCorrection ? "true" : "false") << endl;
}
