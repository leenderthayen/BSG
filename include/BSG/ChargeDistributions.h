#ifndef CHARGEDISTRIBUTIONS
#define CHARGEDISTRIBUTIONS

#include "Utilities.h"
#include "gsl/gsl_sf_laguerre.h"
#include "gsl/gsl_sf_gamma.h"

#include <complex>

#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"

namespace bsg {
/**
 * Namespace containing all kinds of utility functions related to
 * charge distributions, typically those with harmonic oscillator
 * functions
*/
namespace ChargeDistributions {

using std::cout;
using std::endl;

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
  //auto start = std::chrono::steady_clock::now();
  TF1* funcMG = new TF1("ChargeMG", ChargeMG_f, 0, 5 * rms, 2);
  funcMG->SetParameters(5.0, std::sqrt(5. / 3.) * rms);
  funcMG->FixParameter(1, std::sqrt(5. / 3.) * rms);

  /*auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "Functions milliseconds since start: " << elapsed.count() << "\n";*/
  Double_t x[50], y[50];
  Int_t n = 50;
  for (Int_t i = 0; i < n; i++) {
    x[i] = i * 5 * rms / n;
    y[i] = ChargeHO(x[i], rms, Z, true);
  }

  TGraph* gr = new TGraph(n, x, y);

  /*elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "Filling milliseconds since start: " << elapsed.count() << "\n";*/

  gr->Fit(funcMG, "WQ0");
  /*elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  std::cout << "Fitting milliseconds since start: " << elapsed.count() << "\n";*/

  double A = funcMG->GetParameter(0);

  delete funcMG;
  //delete funcHO;
  //delete histHO;
  return A;
}
}
}
#endif
