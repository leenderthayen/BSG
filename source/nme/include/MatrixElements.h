#ifndef MATRIXELEMENTS
#define MATRIXELEMENTS

#include "gsl/gsl_sf_coupling.h"
#include "NilssonOrbits.h"
#include "Utilities.h"
#include "ChargeDistributions.h"
#include "NuclearUtilities.h"

#include <cmath>
#include <complex>

namespace NuclearStructure {
/**
 * Namespace collecting all kinds of functiosn related to single particle matrix
 * elements and angular momentum coupling
 */
namespace MatrixElements {

using std::cout;
using std::endl;

namespace CD = ChargeDistributions;

/**
 * Get the orbital angular momentum value from the kappa value
 *
 * @param k kappa value of state
 */
inline int gL(int k) { return (k > 0) ? k : std::abs(k) - 1; }

/**
 * Get the kappa value for an ls state
 *
 * @param l orbital angular momentum
 * @param s spin direction
 */
inline int kL(int l, int s) { return (s > 0) ? -(l + 1) : l; }

/**
 * Get the total angular momentum from the kappa value
 *
 * @param k kappa
 */
inline double jL(int k) { return (k < 0) ? -k - 0.5 : k - 0.5; }

/**
 * Calculate the angular momentum function @f[ G_{KLs}(\kappa_i, \kappa_f) @f]
 * as first defined by Weidenmueller
 *
 * @param ki kappa of initial state
 * @param kf kappa of final state
 * @param K spherical tensor rank of the operator
 * @param L orbital angular momentum of the operator
 * @param s index denoting simple of vector spherical tensor
 */
inline double CalculateGKLs(int kf, int ki, int K, int L, int s) {
  int dJi = 2 * jL(ki);
  int dJf = 2 * jL(kf);
  double first = std::sqrt((2 * s + 1) * (2 * K + 1) * (2 * gL(kf) + 1) *
                           (2 * gL(ki) + 1) * (dJf + 1) * (dJi + 1));
  std::complex<double> I(0, 1);
  std::complex<double> sPow =
      std::pow(I, gL(ki) + gL(kf) + L) * std::pow(-1, (dJi - dJf) / 2.);
  double second = sPow.real();
  double third =
      utilities::ClebschGordan(2 * gL(kf), 2 * gL(ki), 2 * L, 0, 0, 0);
  double fourth = gsl_sf_coupling_9j(2 * K, 2 * s, 2 * L, dJf, 1, 2 * gL(kf),
                                     dJi, 1, 2 * gL(ki));

  return first * second * third * fourth;
}

/**
 * Calculate the spin-reduced single particle matrix element @f[ ^{V/A}M_{KLs} @f]
 * as calculated between two spherical harmonic oscillator states
 * in the non-relativistic approximation. We work in the qR << 1 limit, where
 * @f[ j_L(qr) \approx \frac{(qr)^L}{(2L+1)!!} @f], i.e. the matrix element
 * is insensitive to the momentum transfer
 *
 * @param V boolean denoting vector (true) or axial vector (false) nature of
 *matrix element
 * @param Ji spin of initial state, used for normalization
 * @param K spherical tensor rank of the operator
 * @param L orbital angular momentum of the operator
 * @param s index denoting simple or vector spherical tensor
 * @param ni principal quantum number of initial state
 * @param nf principal quantum number of final state
 * @param li orbital angular momentum of initial state
 * @param lf orbital angular momentum of final state
 * @param si spin direction of initial state
 * @param sf spin direction of final state
 * @param R nuclear radius
 * @param nu length scale of the harmonic oscillator functions
 */
inline double GetReducedSingleParticleMatrixElement(bool V, double Ji, int K, int L,
                                             int s, int ni, int nf, int li,
                                             int lf, int si, int sf, double R,
                                             double nu) {
  double Mn = NUCLEON_MASS_KEV / ELECTRON_MASS_KEV;
  double result = std::sqrt(2. / (2. * Ji + 1.));

  int kf = kL(lf, sf);
  int ki = kL(li, si);

  if (V) {
    if (s == 0) {
      result *= CalculateGKLs(kf, ki, K, L, 0.);
      result *= CD::GetRadialMEHO(nf, lf, K, ni, li, nu) / std::pow(R, K);
    } else if (s == 1) {
      double dE = 2. * nu * ELECTRON_MASS_KEV / NUCLEON_MASS_KEV *
                  (2 * (ni - nf) + li - lf);
      double first =
          R / 2. / (L + 1.) * dE *
              CD::GetRadialMEHO(nf, lf, L + 1, ni, li, nu) /
              std::pow(R, L + 1) +
          (kf - ki + 1 + L) * (kf + ki - L) /
              (4. * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV * R) / (L + 1.) *
              CD::GetRadialMEHO(nf, lf, L - 1, ni, li, nu) / std::pow(R, L - 1);
      double second =
          -R / 2. / (L + 1.) * dE *
              CD::GetRadialMEHO(nf, lf, L + 1, ni, li, nu) /
              std::pow(R, L + 1) -
          (kf - ki - 1 - L) * (kf + ki - L) /
              (4. * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV * R) / (L + 1.) *
              CD::GetRadialMEHO(nf, lf, L - 1, ni, li, nu) / std::pow(R, L - 1);
      result *= sign(ki) * CalculateGKLs(kf, -ki, K, L, s) * first +
                sign(kf) * CalculateGKLs(-kf, ki, K, L, s) * second;
    } else {
      result = 0.0;
    }
  } else {
    if (s == 0) {
      double dE = 2. * nu * ELECTRON_MASS_KEV / NUCLEON_MASS_KEV *
                  (2 * (ni - nf) + li - lf);
      double first =
          R / 2. / (K + 1.) * dE *
              CD::GetRadialMEHO(nf, lf, K + 1, ni, li, nu) /
              std::pow(R, K + 1) +
          (kf - ki + 1 + K) * (kf + ki - K) /
              (4. * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV * R) / (K + 1.) *
              CD::GetRadialMEHO(nf, lf, K - 1, ni, li, nu) / std::pow(R, K - 1);
      double second =
          -R / 2. / (K + 1.) * dE *
              CD::GetRadialMEHO(nf, lf, K + 1, ni, li, nu) /
              std::pow(R, K + 1) -
          (kf - ki - 1 - K) * (kf + ki - K) /
              (4. * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV * R) / (K + 1.) *
              CD::GetRadialMEHO(nf, lf, K - 1, ni, li, nu) / std::pow(R, K - 1);
      result *= sign(ki) * CalculateGKLs(kf, -ki, K, L, 0) * first +
                sign(kf) * CalculateGKLs(-kf, ki, K, L, s) * second;
    } else if (s == 1) {
      result *= CalculateGKLs(kf, ki, K, L, s);
      result *= CD::GetRadialMEHO(nf, lf, L, ni, li, nu) / std::pow(R, L);
    } else {
      result = 0.0;
    }
  }

  return result;
}

/**
 * Calculate spin-reduced single particle matrix using initial and final single particle
 *states
 * containing a mixture of harmonic oscillator wave functions.
 *
 * @param V boolean denoting vector (true) or axial vector (false) nature of
 *matrix element
 * @param Ji nuclear spin of initial state
 * @param K rank of spherical tensor of the operator
 * @param L orbital angular momentum of the operator
 * @param s index denoting simple or vector spherical harmonics
 * @param spsi initial single particle state
 * @param spsf final single particle state
 * @param R nuclear radius
 * @param nu length scale of the harmonic oscillator functions
 */
inline double GetReducedSingleParticleMatrixElement(bool V, double Ji, int K, int L,
                                             int s, SingleParticleState spsi,
                                             SingleParticleState spsf, double R,
                                             double nu) {
  std::vector<WFComp> initComps = spsi.componentsHO;
  std::vector<WFComp> finalComps = spsf.componentsHO;

  double result = 0.0;

  for (int i = 0; i < initComps.size(); i++) {
    for (int j = 0; j < finalComps.size(); j++) {
      result += initComps[i].C * finalComps[j].C *
                GetReducedSingleParticleMatrixElement(V, Ji, K, L, s, initComps[i].n,
                                               finalComps[j].n, initComps[i].l,
                                               finalComps[j].l, initComps[i].s,
                                               finalComps[j].s, R, nu);
    }
  }

  return result;
}

/**
 * Calculate the phase when switching the value of Omega for the C coefficients
 *
 * @param w wave function component
 * @param dO double of Omega value
 */
inline double GetCjO(WFComp w, int dO) {
  if (dO > 0) {
    return w.C;
  } else {
    return std::pow(-1, w.s / 2. - 0.5) * w.C;
  }
}

/**
 * Calculate the reduced single particle matrix element in a deformed nucleus
 * between two single particle states
 *
 * @param opt option parameter for the difference between odd-A (0), odd-odd ->
 *even-even (1)
 *     and even-even -> odd-odd (2) transition
 * @param spsi initial single particle state
 * @param spsf final single particle state
 * @param V boolean denoting vector (true) or axial vector (false) nature of
 *matrix element
 * @param K rank of spherical tensor of the operator
 * @param L orbital angular momentum of the operator
 * @param s index denoting simple or vector spherical harmonics
 * @param dJi double of initial nuclear spin
 * @param dJf double of final nuclear spin
 * @param dKi double of initial K value
 * @param dKf double of final K value
 * @param R nuclear radius in natural units
 * @param nu length scale of harmonic oscillator wave functions in natural units
 * @see GetSingleParticleMatrixElement
 */
inline double GetDeformedReducedSingleParticleMatrixElement(
    int opt, SingleParticleState spsi, SingleParticleState spsf, bool V, int K,
    int L, int s, int dJi, int dJf, int dKi, int dKf, double R, double nu) {
  double result = 0.;

  // Odd-A transition
  if (opt == 0) {
    std::vector<WFComp> finalStates = spsf.componentsHO;
    std::vector<WFComp> initStates = spsi.componentsHO;
    int inO = spsi.dO;
    int fO = spsf.dO;
    for (int i = 0; i < finalStates.size(); i++) {
      WFComp fW = finalStates[i];
      for (int j = 0; j < initStates.size(); j++) {
        WFComp iW = initStates[j];
        result +=
            fW.C * iW.C *
            (std::pow(-1., (dJf - dKf + 2 * fW.l + fW.s - fO) / 2.) *
                 gsl_sf_coupling_3j(dJf, 2 * K, dJi, -dKf, fO - inO, dKi) *
                 gsl_sf_coupling_3j(2 * fW.l + fW.s, 2 * K, 2 * iW.l + iW.s,
                                    -fO, fO - inO, inO) +
             gsl_sf_coupling_3j(dJf, 2 * K, dJi, dKf, -fO - inO, dKi) *
                 gsl_sf_coupling_3j(2 * fW.l + fW.s, 2 * K, 2 * iW.l + iW.s, fO,
                                    -fO - inO, inO)) *
            GetReducedSingleParticleMatrixElement(V, dJi / 2., K, L, s, iW.n, fW.n,
                                           iW.l, fW.l, iW.s, fW.s, R, nu);
      }
    }

    result *= std::sqrt((dJi + 1.) * (dJf + 1.) / (1. + delta(dKf, 0.)) /
                        (1. + delta(dKi, 0.)));
  } else {
    // Even-A transition
    std::vector<WFComp> finalStates = spsf.componentsHO;
    std::vector<WFComp> initStates = spsi.componentsHO;
    int inO = spsi.dO;
    int fO = spsf.dO;

    if (opt == 1) {
      // Odd-Odd ---> Even-Even
      for (int i = 0; i < finalStates.size(); i++) {
        WFComp fW = finalStates[i];
        for (int j = 0; j < initStates.size(); j++) {
          WFComp iW = initStates[j];
          result +=
              GetCjO(fW, -fO) * GetCjO(iW, inO) *
              std::pow(-1., fW.l + fW.s / 2. + fO / 2.) *
              gsl_sf_coupling_3j(2 * fW.l + fW.s, 2 * K, 2 * iW.l + iW.s, fO,
                                 -dKi, inO) *
              GetReducedSingleParticleMatrixElement(V, dJi / 2., K, L, s, iW.n, fW.n,
                                             iW.l, fW.l, iW.s, fW.s, R, nu);
        }
      }
      // Even-Even ---> Odd-Odd
    } else if (opt == 2) {
      cout << fO << " " << dKf << " " << inO << endl;
      for (int i = 0; i < finalStates.size(); i++) {
        WFComp fW = finalStates[i];
        for (int j = 0; j < initStates.size(); j++) {
          WFComp iW = initStates[j];
          result +=
              GetCjO(fW, fO) * GetCjO(iW, -inO) *
              std::pow(-1., fW.l + fW.s / 2. - fO / 2.) *
              gsl_sf_coupling_3j(2 * fW.l + fW.s, 2 * K, 2 * iW.l + iW.s, -fO,
                                 dKf, -inO) *
              GetReducedSingleParticleMatrixElement(V, dJi / 2., K, L, s, iW.n, fW.n,
                                             iW.l, fW.l, iW.s, fW.s, R, nu);
        }
      }
    }
  }
  return result;
}
}
}
#endif
