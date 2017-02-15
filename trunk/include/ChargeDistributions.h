namespace chargedistributions {
inline double GetRMSHO(int n, int l, double nu) {
  return std::sqrt(1. / 4. / nu * (4 * n + 2 * l - 1));
}

inline double RadialHO(int n, int l, double nu, double r) {
  const int k = n - 1;

  double N =
      std::sqrt(sqrt(2 * std::pow(nu, 3) / M_PI) * std::pow(2., k + 2 * l + 3) * Factorial(k) *
           std::pow(nu, l) / DoubleFactorial(2 * k + 2 * l + 1));

  double gl = gsl_sf_laguerre_n(k, l + 0.5, 2 * nu * r * r);

  return N * std::pow(r, l + 1) * exp(-nu * r * r) * gl;
}

inline double CalcBoverA(double beta) {
  double C = 5./3.*std::sqrt(5./M_PI)*beta*(1.+0.16*beta);
  return std::sqrt((1+C)/(1-C/2.));
}

struct RadialParams {double nu1, nu2; int n1, l1, n2, l2; };

inline double RadialHORMS(double r, void * p) {
  struct RadialParams * params = (RadialParams*)p;

  double nu1 = params->nu1;
  double nu2 = params->nu2;
  int n1 = params->n1;
  int n2 = params->n2;
  int l1 = params->l1;
  int l2 = params->l2;

  //cout << "nu1: " << nu1 << " n2: " << n2 << " l1: " << l1 << endl;

  return r*r*abs(RadialHO(n1, l1, nu1, r)*RadialHO(n2, l2, nu2, r));
}

inline double RadialHONorm(double r, void * p) {
  struct RadialParams * params = (RadialParams*)p;

  double nu1 = params->nu1;
  double nu2 = params->nu2;
  int n1 = params->n1;
  int n2 = params->n2;
  int l1 = params->l1;
  int l2 = params->l2;

  return abs(RadialHO(n1, l1, nu1, r)*RadialHO(n2, l2, nu2, r));
}

inline double WeakIntegratedRMS(double nu1, int n1, int l1, double nu2, int n2, int l2) {
  int intervals = 2000;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(intervals);

  double resultr2, errorr2;
  gsl_function Fr2;
  struct RadialParams params = {nu1, nu2, n1, l1, n2, l2};
  
  Fr2.function = &RadialHORMS;
  Fr2.params = &params;

  double epsabs = 1.0E-4;
  double epsrel = 1.0E-4;

  int status = gsl_integration_qag(&Fr2, 0, 0.1, epsabs, epsrel, intervals, 6, w, &resultr2, &errorr2);

  if (status) {
    if (status == GSL_EMAXITER) {
      cout << "ERROR: Maximum number of subdivisions reached." << endl;
    }
    else {
      cout << "ERROR: " << gsl_strerror(status) << endl;
    }
    return 0.;
  }

  double resultNorm, errorNorm;
  gsl_function FNorm;

  FNorm.function = &RadialHONorm;
  FNorm.params = &params;

  status = gsl_integration_qag(&FNorm, 0, 0.1, epsabs, epsrel, intervals, 6, w, &resultNorm, &errorNorm);

  gsl_integration_workspace_free(w);

  return resultr2/resultNorm;
}

inline double CalcNu(double rms, int Z) {
  std::vector<int> occNumbers = GetOccupationNumbers(Z);

  double s = 0.;
  for (int i = 0; i < occNumbers.size(); i += 4) {
    s += (4 * occNumbers[i] + 2 * occNumbers[i + 1] - 1) * occNumbers[i+3];
  }
  return s / Z / 4. / rms / rms;
}

inline double CalcChargeIndepNu(double rms, int Z, int N) {
  std::vector<int> occNumbersZ = GetOccupationNumbers(Z);
  std::vector<int> occNumbersN = GetOccupationNumbers(N);

  double s = 0.;
  for (int i = 0; i < occNumbersZ.size(); i += 4) {
    s += (4 * occNumbersZ[i] + 2 * occNumbersZ[i + 1] - 1) * occNumbersZ[i+3];
  }
  for (int i = 0; i < occNumbersN.size(); i += 4) {
    s += (4 * occNumbersN[i] + 2 * occNumbersN[i + 1] - 1) * occNumbersN[i+3];
  }
  return s / (Z + N) / 4. / rms / rms;
}

inline double ChargeHO(double r, double rms, int Z, bool normalised) {
  double nu = CalcNu(rms, Z);

  double charge = 0.;

  std::vector<int> occNumbers = GetOccupationNumbers(Z);
  for (int i = 0; i < occNumbers.size(); i += 4) {
    charge += std::pow(RadialHO(occNumbers[i], occNumbers[i + 1], nu, r), 2.) *
              occNumbers[i + 3];
  }
  if (normalised) {
    charge /= Z;
  }

  return charge;
}

inline double ChargeMG(double r, double A, double R) {
  double a = R/std::sqrt(5./2.*(2.+5.*A)/(2.+3.*A));
  return 8./(2.+3.*A)/std::pow(a, 3)/std::sqrt(M_PI)*(1.+A*std::pow(r/a, 2.))*exp(-std::pow(r/a, 2.))*r*r;
}

inline Double_t ChargeHO_f(Double_t x[], Double_t par[]) {
  Double_t f;

  f = ChargeHO(x[0], par[0], (int)par[1], true);

  return f;
}

inline Double_t ChargeMG_f(Double_t x[], Double_t par[]) {
  Double_t f;

  f = ChargeMG(x[0], par[0], par[1]);

  return f;
}

inline double DeformedChargeDist(double r, double a, double b) {
  double rho0 = 3./4./M_PI/a/a/b;
  double result = 0.;
  if (r <= a) {
    result = rho0;
  }
  else if (r <= b) {
    result = rho0*(1-b/r*std::sqrt((r*r-a*a)/(b*b-a*a)));
  }
  return result;
}

inline double GetDerivDeformedChargeDist(double r, double a, double b) {
  double rho0 = 3./4./M_PI/a/a/b;
  double result = 0.;
  if (r >= a && r <= b) {
    result = rho0*(b*std::sqrt((r*r-a*a)/(b*b-a*a))/r/r-b/std::sqrt((b*b-a*a)*(r*r-a*a)));
  }
  return result;
}

inline double FitHODist(int Z, double rms) {

  TF1* funcMG = new TF1("ChargeMG", ChargeMG_f, 0, 5*rms, 2);
  funcMG->SetParameters(2.0, std::sqrt(5./3.)*rms);
  funcMG->FixParameter(1, std::sqrt(5./3.)*rms);

  TF1* funcHO = new TF1("ChargeHO", ChargeHO_f, 0, 5*rms, 2);
  funcHO->SetParameters(rms, (double)Z);
  TH1F* histHO = new TH1F("stats", "my pdf", 1000, 0, 5*rms);

  double N = 1E7;
  for(Int_t i = 0; i<N; i++) {
    histHO->Fill(funcHO->GetRandom());
  }
  histHO->Scale(1./histHO->Integral(), "width");
  histHO->Fit("ChargeMG", "WQ0");

  double A = funcMG->GetParameter(0);

  delete funcMG;
  delete funcHO;
  delete histHO;
  return A;
}
}
