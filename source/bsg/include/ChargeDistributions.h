#ifndef CHARGEDISTRIBUTIONS
#define CHARGEDISTRIBUTIONS

#include "Utilities.h"
#include "gsl/gsl_sf_laguerre.h"
#include "gsl/gsl_sf_gamma.h"

#include <complex>

#include "TMath.h"
#include "TF1.h"
#include "TH1.h"

/**
 * Namespace containing all kinds of utility functions related to 
 * charge distributions, typically those with harmonic oscillator
 * functions 
*/
namespace ChargeDistributions {

using std::cout;
using std::endl;

/**
 * Calculates the root mean square (RMS) radius of
 * a harmonic oscillator (HO) function defined by nu.
 *
 * @param n Main radial quantum number
 * @param l orbital quantum number
 * @param nu Length scale of the HO function
 * @returns the RMS radius in the legnth scale of nu
 * @see GetRadialMEHO
 */
inline double GetRMSHO(int n, int l, double nu) {
  return std::sqrt(1. / 4. / nu * (4 * n + 2 * l - 1));
}

/**
 * Get the general radial matrix element for a harmonic oscillator function.
 * @f[\langle n_fl_l | r^L | n_il_i \rangle @f]
 * 
 * @param nf the main radial quantum number of the final state
 * @param lf the orbital quantum number of the final state
 * @param exponent of the radial matrix element
 * @param ni the main radial quantum number of the initial state
 * @param li the orbital quantum number of the initial state
 * @param nu the length scale of the harmonic oscillator function
 * @returns the radial main element of power L between two HO states
 */
inline double GetRadialMEHO(int nf, int lf, int L, int ni, int li, double nu) {
  int taui = (lf - li + L) / 2;
  int tauf = (li - lf + L) / 2;
  double t = (li + lf + L + 1) / 2.;

  double first =
      std::pow(-1, ni + nf) / std::pow(2 * nu, L / 2.) *
      std::sqrt(gsl_sf_gamma(ni) * gsl_sf_gamma(nf) /
                gsl_sf_gamma(ni + t - taui) / gsl_sf_gamma(nf + t - tauf)) *
      utilities::Factorial(taui) * utilities::Factorial(tauf);

  double sum = 0.;
  for (int s = std::max(std::max(ni - taui - 1, nf - tauf - 1), 0);
       s <= std::min(ni - 1, nf - 1); s++) {
    sum += gsl_sf_gamma(t + s + 1.) /
           (utilities::Factorial(s) * utilities::Factorial(ni - s - 1) *
            utilities::Factorial(s + taui - ni + 1) *
            utilities::Factorial(s + tauf - nf + 1));
  }

  return first * sum;
}

/**
 * Gives the value of a harmonic oscillator function at radius r,
 * defined as the generalized Laguerre polynomial
 * 
 * @param n the main radial quantum number of the final state
 * @param l the orbital quantum number
 * @param nu the length scale of the harmonic oscillator function
 * @param r the radius, in the same length scale as nu
 * @returns the value of the HO function at radius r
 */
inline double RadialHO(int n, int l, double nu, double r) {
  const int k = n - 1;

  double N =
      std::sqrt(sqrt(2 * std::pow(nu, 3) / M_PI) * std::pow(2., k + 2 * l + 3) *
                utilities::Factorial(k) * std::pow(nu, l) /
                utilities::DoubleFactorial(2 * k + 2 * l + 1));

  double gl = gsl_sf_laguerre_n(k, l + 0.5, 2 * nu * r * r);

  return N * std::pow(r, l + 1) * exp(-nu * r * r) * gl;
}

/**
 * Calculate the ratio of long over short edge of an ellipse
 * based on the nuclear quadrupole deformation
 * 
 * @param beta the nuclear quadrupole deformation parameter
 * @returns the ratio of the long and short edge of the ellipsoid
 */
inline double CalcBoverA(double beta) {
  double C = 5. / 3. * std::sqrt(5. / M_PI) * beta * (1. + 0.16 * beta);
  return std::sqrt((1 + C) / (1 - C / 2.));
}

/**
 * Helper struct in calculating the overlap integral between two
 * HO states.
 */
struct RadialParams {
  double nu1, nu2;
  int n1, l1, n2, l2;
};

/**
 * 
 * 
 * @param r the radius
 * @param p void pointer to a RadialParams struct
 * @returns the value of r^2 weighted with both HO functions
 */
inline double RadialHORMS(double r, void* p) {
  struct RadialParams* params = (RadialParams*)p;

  double nu1 = params->nu1;
  double nu2 = params->nu2;
  int n1 = params->n1;
  int n2 = params->n2;
  int l1 = params->l1;
  int l2 = params->l2;

  // cout << "nu1: " << nu1 << " n2: " << n2 << " l1: " << l1 << endl;

  return r * r * abs(RadialHO(n1, l1, nu1, r) * RadialHO(n2, l2, nu2, r));
}

/**
 * 
 * 
 * @param r the radius
 * @param p void pointer to a RadialParams struct
 * @returns the absolute value of the radial overlap between two HO states
 */
inline double RadialHONorm(double r, void* p) {
  struct RadialParams* params = (RadialParams*)p;

  double nu1 = params->nu1;
  double nu2 = params->nu2;
  int n1 = params->n1;
  int n2 = params->n2;
  int l1 = params->l1;
  int l2 = params->l2;

  return abs(RadialHO(n1, l1, nu1, r) * RadialHO(n2, l2, nu2, r));
}

/** 
 * Calculated the `weak' RMS value of two harmonic oscillator states.
 * 
 * @param nu1 length scale of the first HO function
 * @param n1 principal quantum number of the first HO state
 * @param l1 orbital quantum number of the first HO state
 * @param nu2 length scale of the second HO function
 * @param n2 principal quantum number of the second HO state
 * @param l2 orbital quantum number of the second HO state
 * @returns the normalized RMS radius of two HO states through explicit integration
 */
inline double WeakIntegratedRMS(double nu1, int n1, int l1, double nu2, int n2,
                                int l2) {
  int intervals = 2000;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc(intervals);

  double resultr2, errorr2;
  gsl_function Fr2;
  struct RadialParams params = {nu1, nu2, n1, l1, n2, l2};

  Fr2.function = &RadialHORMS;
  Fr2.params = &params;

  double epsabs = 1.0E-4;
  double epsrel = 1.0E-4;

  int status = gsl_integration_qag(&Fr2, 0, 0.1, epsabs, epsrel, intervals, 6,
                                   w, &resultr2, &errorr2);

  if (status) {
    if (status == GSL_EMAXITER) {
      cout << "ERROR: Maximum number of subdivisions reached." << endl;
    } else {
      cout << "ERROR: " << gsl_strerror(status) << endl;
    }
    return 0.;
  }

  double resultNorm, errorNorm;
  gsl_function FNorm;

  FNorm.function = &RadialHONorm;
  FNorm.params = &params;

  status = gsl_integration_qag(&FNorm, 0, 0.1, epsabs, epsrel, intervals, 6, w,
                               &resultNorm, &errorNorm);

  gsl_integration_workspace_free(w);

  return resultr2 / resultNorm;
}

/**
 * Calculate the length scale nu of a HO state using the nuclear RMS radius and its Z value
 * by looking at the occupation numbers in a spherical harmonic oscillator way through jj-coupling
 * 
 * @param rms the RMS nuclear radius
 * @param Z the proton number of the nucleus
 * @return the appropriate nu value 
 */
inline double CalcNu(double rms, int Z) {
  std::vector<int> occNumbers = utilities::GetOccupationNumbers(Z);

  double s = 0.;
  for (int i = 0; i < occNumbers.size(); i += 4) {
    s += (4 * occNumbers[i] + 2 * occNumbers[i + 1] - 1) * occNumbers[i + 3];
  }
  return s / Z / 4. / rms / rms;
}

/**
 * Calculate the length scale nu of the HO states assuming neutrons and protons both to contribute.
 * 
 * @param rms the RMS nuclear radius
 * @param Z the proton number of the nucleus
 * @returns the appropriate nu value
 * @see CalcNu
 */
inline double CalcChargeIndepNu(double rms, int Z, int N) {
  std::vector<int> occNumbersZ = utilities::GetOccupationNumbers(Z);
  std::vector<int> occNumbersN = utilities::GetOccupationNumbers(N);

  double s = 0.;
  for (int i = 0; i < occNumbersZ.size(); i += 4) {
    s += (4 * occNumbersZ[i] + 2 * occNumbersZ[i + 1] - 1) * occNumbersZ[i + 3];
  }
  for (int i = 0; i < occNumbersN.size(); i += 4) {
    s += (4 * occNumbersN[i] + 2 * occNumbersN[i + 1] - 1) * occNumbersN[i + 3];
  }
  return s / (Z + N) / 4. / rms / rms;
}

/**
 * Calculate the nuclear charge density based on harmonic oscillator states
 * 
 * @param r the radius
 * @param rms the nuclear RMS radius
 * @param Z the proton number of the nucleus
 * @param normalised whether to normalize to unit charge
 * @returns the charge density at radius r
 */
inline double ChargeHO(double r, double rms, int Z, bool normalised) {
  double nu = CalcNu(rms, Z);

  double charge = 0.;

  std::vector<int> occNumbers = utilities::GetOccupationNumbers(Z);
  for (int i = 0; i < occNumbers.size(); i += 4) {
    charge += std::pow(RadialHO(occNumbers[i], occNumbers[i + 1], nu, r), 2.) *
              occNumbers[i + 3];
  }
  if (normalised) {
    charge /= Z;
  }

  return charge;
}

/**
 * Calculate the nuclear charge density based on a modified gaussian distribution
 * 
 * @param r the radius
 * @param A the parameter of the modified gaussian distribution
 * @param R the nuclear radius
 * @returns the charge density at radius r
 */
inline double ChargeMG(double r, double A, double R) {
  double a = R / std::sqrt(5. / 2. * (2. + 5. * A) / (2. + 3. * A));
  return 8. / (2. + 3. * A) / std::pow(a, 3) / std::sqrt(M_PI) *
         (1. + A * std::pow(r / a, 2.)) * exp(-std::pow(r / a, 2.)) * r * r;
}

/**
 * Wrapper function around the Harmonic Oscillator charge density functions used by the ROOT fitter
 * 
 * @param x array of radii
 * @param par vector object containing all other information required by the function
 * @returns the nuclear charge density at radius x[0]
 */
inline Double_t ChargeHO_f(Double_t x[], Double_t par[]) {
  Double_t f;

  f = ChargeHO(x[0], par[0], (int)par[1], true);

  return f;
}

/**
 * Wrapper function around the Modified Gaussian charge density functions used by the ROOT fitter
 * 
 * @param x array of radii
 * @param par vector object containing all other information required by the function
 * @returns the nuclear charge density at x[0]
 */
inline Double_t ChargeMG_f(Double_t x[], Double_t par[]) {
  Double_t f;

  f = ChargeMG(x[0], par[0], par[1]);

  return f;
}

/**
 * Calculate the angle-averaged charge density distribution for a quadrupole deformed nucleus
 * 
 * @param r the radius
 * @param a the short radius in the ellipsoid
 * @param b the long radius in the ellipsoid
 * @returns the charge density at radius r
 */
inline double DeformedChargeDist(double r, double a, double b) {
  double rho0 = 3. / 4. / M_PI / a / a / b;
  double result = 0.;
  if (r <= a) {
    result = rho0;
  } else if (r <= b) {
    result = rho0 * (1 - b / r * std::sqrt((r * r - a * a) / (b * b - a * a)));
  }
  return result;
}

/**
 * Calculate the derivative of the angle-averaged charge density distribution for a quadrupole
 * deformed nucleus
 * 
 * @param r the radius
 * @param a the short radius in the ellipsoid
 * @param b the long radius in the ellipsoid
 * @returns the derivative of the charge density at radius r
 * @see DeformedChargeDist
 */
inline double GetDerivDeformedChargeDist(double r, double a, double b) {
  double rho0 = 3. / 4. / M_PI / a / a / b;
  double result = 0.;
  if (r >= a && r <= b) {
    result = rho0 * (b * std::sqrt((r * r - a * a) / (b * b - a * a)) / r / r -
                     b / std::sqrt((b * b - a * a) * (r * r - a * a)));
  }
  return result;
}

/**
 * Fit the charge distribution as constructed using Harmonic Oscillator functions
 * to a Modified Gaussian distribution using ROOT.
 * 
 * @param Z the proton number of the nucleus
 * @param rms the nuclear RMS radius
 * @returns the fitted A parameter of the Modified Gaussian density
 * @see ChargeHO_f
 * @see ChargeMG_f
 */
inline double FitHODist(int Z, double rms) {
  TF1* funcMG = new TF1("ChargeMG", ChargeMG_f, 0, 5 * rms, 2);
  funcMG->SetParameters(2.0, std::sqrt(5. / 3.) * rms);
  funcMG->FixParameter(1, std::sqrt(5. / 3.) * rms);

  TF1* funcHO = new TF1("ChargeHO", ChargeHO_f, 0, 5 * rms, 2);
  funcHO->SetParameters(rms, (double)Z);
  TH1F* histHO = new TH1F("stats", "my pdf", 1000, 0, 5 * rms);

  double N = 1E7;
  for (Int_t i = 0; i < N; i++) {
    histHO->Fill(funcHO->GetRandom());
  }
  histHO->Scale(1. / histHO->Integral(), "width");
  histHO->Fit("ChargeMG", "WQ0");

  double A = funcMG->GetParameter(0);

  delete funcMG;
  delete funcHO;
  delete histHO;
  return A;
}
}
#endif
