#ifndef MATRIXELEMENTS
#define MATRIXELEMENTS

#include "gsl/gsl_sf_coupling.h"
#include "NilssonOrbits.h"
#include "Utilities.h"
#include "ChargeDistributions.h"

#include <cmath>
#include <complex>

namespace MatrixElements {

namespace CD = ChargeDistributions;

inline int gL(int k) { return (k > 0) ? k : std::abs(k) - 1; }

inline int kL(int l, int s) { return (s > 0) ? -(l + 1) : l; }

inline double jL(int k) { return (k > 0) ? -k - 0.5 : k - 0.5; }

inline double sign(double x) { return (x == 0) ? 0 : x / std::abs(x); }

inline double CalculateGKLs(int kf, int ki, int K, int L, int s) {
  int dJi = 2 * jL(ki);
  int dJf = 2 * jL(kf);
  // cout << "Calc GKLs " << kf << " " << ki << endl;
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

inline double GetSingleParticleMatrixElement(bool V, double Ji, int K, int L,
                                             int s, int ni, int nf, int li,
                                             int lf, int si, int sf, double R,
                                             double nu) {
  // cout << "Calculating SP ME " << K << " " << L << " " << s << " state: " <<
  // li << " " << si << " " << lf << " " << sf << endl;
  double Mn = nucleonMasskeV / electronMasskeV;
  double result = std::sqrt(2. / (2. * Ji + 1.));

  int kf = kL(lf, sf);
  int ki = kL(li, sf);

  if (V) {
    if (s == 0) {
      result *= CalculateGKLs(kf, ki, K, L, 0.);
      result *= CD::GetRadialMEHO(nf, lf, K, ni, li, nu) / std::pow(R, K);
    } else if (s == 1) {
      //TODO check if right energy scale for nu
      double dE = 2. * nu * electronMasskeV / nucleonMasskeV *
                  (2 * (ni - nf) + li - lf);
      double first =
          R / 2. / (L + 1.) * dE * CD::GetRadialMEHO(nf, lf, L + 1, ni, li, nu) /
              std::pow(R, L + 1) +
          (kf - ki + 1 + L) * (kf + ki - L) /
              (4. * nucleonMasskeV / electronMasskeV * R) / (L + 1.) *
              CD::GetRadialMEHO(nf, lf, L - 1, ni, li, nu) / std::pow(R, L - 1);
      double second =
          -R / 2. / (L + 1.) * dE * CD::GetRadialMEHO(nf, lf, L + 1, ni, li, nu) /
              std::pow(R, L + 1) -
          (kf - ki - 1 - L) * (kf + ki - L) /
              (4. * nucleonMasskeV / electronMasskeV * R) / (L + 1.) *
              CD::GetRadialMEHO(nf, lf, L - 1, ni, li, nu) / std::pow(R, L - 1);
      result *= sign(ki) * CalculateGKLs(kf, -ki, K, L, s) * first +
                sign(kf) * CalculateGKLs(-kf, ki, K, L, s) * second;
    } else {
      result = 0.0;
    }
  } else {
    if (s == 0) {
      double dE = 2. * nu * electronMasskeV / nucleonMasskeV *
                  (2 * (ni - nf) + li - lf);
      double first =
          R / 2. / (K + 1.) * dE * CD::GetRadialMEHO(nf, lf, K + 1, ni, li, nu) /
              std::pow(R, K + 1) +
          (kf - ki + 1 + K) * (kf + ki - K) /
              (4. * nucleonMasskeV / electronMasskeV * R) / (K + 1.) *
              CD::GetRadialMEHO(nf, lf, K - 1, ni, li, nu) / std::pow(R, K - 1);
      double second =
          -R / 2. / (K + 1.) * dE * CD::GetRadialMEHO(nf, lf, K + 1, ni, li, nu) /
              std::pow(R, K + 1) -
          (kf - ki - 1 - K) * (kf + ki - K) /
              (4. * nucleonMasskeV / electronMasskeV * R) / (K + 1.) *
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

  if (V) {
    if (s == 0) {
      //result *= getGKLs(kL(lf, sf), kL(li, si), Ji, K, L, s, li, lf, si, sf);
    } else if (s == 1) {
      if (li == lf) {
        if (si == sf) {
          if (si == 1) {
            result *= -(li + 1) * std::sqrt(6. * (li + 1) * (2. * li + 3) /
                                            (2. * li + 1.)) /
                      (2. * Mn * R);
          } else {
            result *= -li *
                      std::sqrt(6. * li * (2 * li - 1.) / (2. * li + 1.)) /
                      (2. * Mn * R);
          }
        } else {
          result *= sf * std::sqrt(3. * li * (li + 1.) / 2. / (2. * li + 1)) /
                    (Mn * R);
        }
      } else {
        result = 0.;
      }
      /*result *=
      1./((protonMasskeV+neutronMasskeV)/electronMasskeV*R)*(sign(kL(li,
      si))*getGKLs(kL(lf, sf), -kL(li, si), Ji, K, L, s, li, lf, si, sf)
                + sign(kL(lf, sf))*getGKLs(-kL(lf, sf), kL(li, si), Ji, K, L, s,
      li, lf, si, sf));
      result *= 3./4.;*/
    }
  } else {
    // result *= getGKLs(kL(lf, sf), kL(li, si), Ji, K, L, s, li, lf, si, sf);
    if (li == lf) {
      if (L == 0) {
        if (si == sf) {
          if (si == 1) {
            result *= std::sqrt((li + 1.) * (2. * li + 3.) / (2. * li + 1));
          } else {
            result *= -std::sqrt(li * (2. * li - 1.) / (2. * li + 1.));
          }
        } else {
          result *= si * 2. * std::sqrt(li * (li + 1.) / (2. * li + 1.));
        }
      } else if (L == 2) {
        if (si == sf) {
          if (si == 1) {
            result *= -li * std::sqrt(2. * (li + 1.) / (2. * li + 1) /
                                      (2. * li + 3.));
          } else {
            result *=
                (li + 1.) * std::sqrt(li / (2. * li + 1) / (2. * li - 1.));
          }
        } else {
          result *= si * std::sqrt(li * (li + 1.) / 2. / (2. * li + 1.));
        }
      }
    } else {
      result = 0.;
    }
  }
  // cout << "Result: " << result << endl;
  return result;
}

inline int delta(double x, double y) { return (x == y) ? 1 : 0; }

inline double CalculateDeformedSPMatrixElement(
    nilsson::SingleParticleState spsi, nilsson::SingleParticleState spsf,
    bool V, int K, int L, int s, double Ji, double Jf, double Ki, double Kf,
    double R, double nu) {
  double result = 0.;

  // cout << "Calculating deformed SP ME " << K << " " << L << " " << s << endl;
  std::vector<nilsson::WFComp> finalStates = spsf.componentsHO;
  std::vector<nilsson::WFComp> initStates = spsi.componentsHO;
  int inO = 2 * spsi.O;
  int fO = 2 * spsf.O;
  for (std::vector<nilsson::WFComp>::iterator fIt = finalStates.begin();
       fIt != finalStates.end(); ++fIt) {
    for (std::vector<nilsson::WFComp>::iterator inIt = initStates.begin();
         inIt != initStates.end(); ++inIt) {
      result +=
          (*fIt).C * (*inIt).C *
          (pow(-1, Jf - Kf + (*fIt).l + (*fIt).s / 2. - fO / 2.) *
               gsl_sf_coupling_3j(2 * Jf, 2 * K, 2 * Ji, -2 * Kf, fO - inO,
                                  2 * Ki) *
               gsl_sf_coupling_3j(2 * (*fIt).l + (*fIt).s, 2 * K,
                                  2 * (*inIt).l + (*inIt).s, -fO, fO - inO,
                                  inO) +
           gsl_sf_coupling_3j(2 * Jf, 2 * K, 2 * Ji, 2 * Kf, -fO - inO,
                              2 * Ki) *
               gsl_sf_coupling_3j(2 * (*fIt).l + (*fIt).s, 2 * K,
                                  2 * (*inIt).l + (*inIt).s, fO, -fO - inO,
                                  inO)) *
          GetSingleParticleMatrixElement(V, Ji, K, L, s, (*inIt).n, (*fIt).n, (*inIt).l, (*fIt).l,
                                         (*inIt).s, (*fIt).s, R, nu);
    }
  }

  result *= std::sqrt((2 * Ji + 1) * (2 * Jf + 1) / (1. + delta(Kf, 0.)) /
                      (1. + delta(Ki, 0.)));

  // cout << "Result: " << result << endl;

  return result;
}
}
#endif
