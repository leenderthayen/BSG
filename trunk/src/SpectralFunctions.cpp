#include "SpectralFunctions.h"
#include "Utilities.h"
#include <complex>
#include <stdio.h>

// GNU Scientific Library stuff
// http://www.gnu.org/software/gsl/
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_const_mksa.h"
#include "gsl/gsl_const_num.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_result.h"
#include "gsl/gsl_sf_dilog.h"

double SpectralFunctions::FermiFunction(double W) {
  double p = sqrt(W * W - 1.);
  double first = 2. * (gamma + 1.);
  // the second term will be incorporated in the fifth
  // double second = 1/pow(gsl_sf_gamma(2*gamma+1),2);
  double third = pow(2. * p * R, 2. * (gamma - 1.));
  double fourth = exp(fBetaType * M_PI * alpha * Zd * W / p);

  // the fifth is a bit tricky
  // we use the complex gamma function from GSL
  gsl_sf_result magn;
  gsl_sf_result phase;
  gsl_sf_lngamma_complex_e(gamma, fBetaType * alpha * Zd * W / p, &magn,
                           &phase);
  // now we have what we want in magn.val

  // but we incorporate the second term here as well
  double fifth = exp(2. * (magn.val - gsl_sf_lngamma(2. * gamma + 1.)));

  // return 1;
  double result = first * third * fourth * fifth;
  return result;
}

double SpectralFunctions::CCorrection(double W) {
  double AC0, AC1, ACm1, AC2;

  /*AC0 = -233. * pow(alpha * Zd, 2) / 630. - (W0 * W0 - 1) * R * R / 5. +
        fBetaType * 2. * W0 * R * alpha * Zd / 35.;
  AC1 = -fBetaType * 21. * R * alpha * Zd / 35.0 + 4.0 / 9.0 * W0 * R * R;
  ACm1 =
      -2 / 45. * gamma * W0 * R * R - fBetaType / 70. * gamma * alpha * Zd * R;
  AC2 = -4. / 9. * R * R;*/

  double VC0, VC1, VCm1, VC2;

  /*VC0 = -233 * pow(alpha * Zd, 2.) / 630. - pow(W0 * R, 2) / 5. -
        fBetaType * 6. / 35. * alpha * Zd * W0 * R +
        9./1225.*pow(alpha*Zd*W0*R, 2.);
  VC1 = -fBetaType * 13. / 35. * alpha * Zd * R + 4. / 15. * W0 * R * R;
  VCm1 = 2. / 15. * gamma * W0 * R * R +
         fBetaType / 70. * gamma * alpha * Zd * R;
  VC2 = -4. / 15. * R * R;*/

  double A;
  if (Afit == -1) {
    double A = utilities::FitHODist(Zd, R * sqrt(3. / 5.));
    Afit = A;
  }

  double F1111 = 0.757 + 0.0069 * (1 - exp(-A / 1.008));
  double F1221 = 0.844 - 0.0182 * (1 - exp(-A / 1.1974));
  double F1222 = 1.219 - 0.0640 * (1 - exp(-A / 1.550));
  double F1211 = -3. / 70;

  VC0 = -pow(W0 * R, 2.) / 5. -
        fBetaType * 2. / 9. * alpha * Zd * W0 * R * F1111 -
        pow(alpha * Zd, 2.) / 3. * F1222;

  VC1 = 4. / 15. * W0 * R * R -
        fBetaType * 2. / 3. * alpha * Zd * R * (F1221 - F1111 / 3.);

  VCm1 = 2. / 15. * W0 * R * R - fBetaType * alpha * Zd * R / 3. * F1211;

  VC2 = -4. / 15. * R * R;

  AC0 = -1. / 3. * pow(alpha * Zd, 2.) * F1222 -
        1. / 5. * (W0 * W0 - 1) * R * R +
        fBetaType * 2. / 27. * alpha * Zd * W0 * R * F1111 + 11. / 45. * R * R;

  AC1 = 4. / 9. * W0 * R * R -
        fBetaType * 2. / 3. * alpha * Zd * R * (1. / 9. * F1111 + F1221);

  ACm1 = -2. / 45. * W0 * R * R + fBetaType * alpha * Zd / 3. * F1211;

  AC2 = -4. / 9. * R * R;

  double result = 1.;

  if (fDecayType == FERMI) {
    result = 1. + VC0 + VC1 * W + VCm1 / W + VC2 * W * W;
  } else if (fDecayType == GAMOW_TELLER) {
    result = 1. + AC0 + AC1 * W + ACm1 / W + AC2 * W * W;
  }
  // TODO
  /*} else if (MixingRatio > 0) {
    result = 1. +
           1. / (1 + pow(MixingRatio, 2)) *
               (VC0 + VC1 * W + VCm1 / W + VC2 * W * W) +
           1. / (1 + pow(MixingRatio, -2)) * (AC0 + AC1 * W + ACm1 / W + AC2 * W
  * W);
  }*/

  if (fCICorrection) {
    result *= CICorrection(W);
  }

  if (fDecayType == GAMOW_TELLER) {
    double M = An * (protonMasskeV + neutronMasskeV) / 2. / electronMasskeV;

    // TODO depends on beta+/- to check from parent or daughter nucleus, not
    // both parent
    int nN, lN, sN, nZ, lZ, sZ;
    vector<int> occNumbersN =
        utilities::GetOccupationNumbers(An - (Zd - fBetaType));
    nN = occNumbersN[occNumbersN.size() - 1 - 3];
    lN = occNumbersN[occNumbersN.size() - 1 - 2];
    sN = occNumbersN[occNumbersN.size() - 1 - 1];
    vector<int> occNumbersZ = utilities::GetOccupationNumbers(Zd - fBetaType);
    nZ = occNumbersZ[occNumbersZ.size() - 1 - 3];
    lZ = occNumbersZ[occNumbersZ.size() - 1 - 2];
    sZ = occNumbersZ[occNumbersZ.size() - 1 - 1];

    // TODO replace <r^2>/R^2 with the HO rms
    double ratioM121 = 0.;
    double nu = utilities::CalcNu(sqrt(3. / 5.) * R, Zd - fBetaType);
    double rHO = pow(utilities::GetRMSHO(nN, lN, nu), 2.);
    /*if (lN == lZ) {
      double nu = utilities::CalcNu(sqrt(3./5.)*R, Zd-fBetaType);
      double rHO = pow(utilities::GetRMSHO(nN, lN, nu), 2.);
      if (sN == sZ) {
        if (sN == 1) {
          ratioM121 = - lN*sqrt(2)/(2.*lN+3.)*rHO/R/R;
        }
        else if (sN == -1) {
          ratioM121 = - (lN+1.)*sqrt(2)/(2.*lN-1.)*rHO/R/R;
        }
      }
      else {
        ratioM121 = 1./2./sqrt(2.)*rHO/R/R;
      }
    }*/

    struct utilities::NuclearState nsi =
        utilities::CalculateDeformedState(Zd - fBetaType, An, beta2);
    struct utilities::NuclearState nsf =
        utilities::CalculateDeformedState(Zd, An, beta2);

    ratioM121 = utilities::CalculateDeformedSPMatrixElement(
                    nsi.states, nsf.states, false, 1, 2, 1, nsi.O, nsf.O, nsi.K,
                    nsf.K, R) /
                utilities::CalculateDeformedSPMatrixElement(
                    nsi.states, nsf.states, false, 1, 0, 1, nsi.O, nsf.O, nsi.K,
                    nsf.K, R);

    ratioM121 *= rHO / R / R;

    double x = sqrt(2.) / 3. * 10. * ratioM121 -
               gP / gA * 25. / 3. / pow(1837.5 * R, 2.);

    double NSC0 = -1. / 45. * R * R * x +
                  1. / 3. * W0 / M / fc1 * (-fBetaType * 2. * fb + fd) +
                  fBetaType * 2. / 5. * alpha * Zd / M / R / fc1 *
                      (fBetaType * 2. * fb + fd) -
                  fBetaType * 2. / 35. * alpha * Zd * W0 * R * x;

    double NSC1 = fBetaType * 4. / 3. * fb / M / fc1 -
                  2. / 45. * W0 * R * R * x +
                  fBetaType * alpha * Zd * R * 2. / 35. * x;

    double NSCm1 = -1. / 3. / M / fc1 * (fBetaType * 2. * fb + fd) +
                   2. / 45. * W0 * R * R * x;

    double NSC2 = 2. / 45. * R * R * x;

    result += NSC0 + NSC1 * W + NSCm1 / W + NSC2 * W * W;
  }

  // cout << "Mixing ratio badly defined. Returning 1." << endl;
  return result;
}

/*double SpectralFunctions::QCDInducedCorrection(double W) {
  double M = An * (protonMasskeV + neutronMasskeV) / 2. / electronMasskeV;
  double h1tildec =
      1 - 2 * W0 / 3. / M / fc * (fd + fBetaType * fb) +
      fBetaType * 4 / 3. * W / M / fc * fb -
      1. / 3 / W / M / fc * (fd + fBetaType * 2 * fb - fh * (W0 - W) / 2 / M);

  double dh1c =
      sqrt(10) / 6. * alpha * Zd / M / R / fc * (2 * fb + fBetaType * fd);
  dh1c = 0.;

  return h1tildec + dh1c;
}*/

double SpectralFunctions::CICorrection(double W) {
  double AC0, AC1, AC2, ACm1;
  double VC0, VC1, VC2, VCm1;

  double rms = sqrt(3. / 5.) * R;
  double nu = 0.;

  int nN, lN, nZ, lZ;
  vector<int> occNumbersN =
      utilities::GetOccupationNumbers(An - (Zd - fBetaType));
  nN = occNumbersN[occNumbersN.size() - 1 - 3];
  lN = occNumbersN[occNumbersN.size() - 1 - 2];
  vector<int> occNumbersZ = utilities::GetOccupationNumbers(Zd - fBetaType);
  nZ = occNumbersZ[occNumbersZ.size() - 1 - 3];
  lZ = occNumbersZ[occNumbersZ.size() - 1 - 2];

  // cout << "nZ: " << nZ << " lZ: " << lZ << " p: " <<
  // occNumbersZ[occNumbersZ.size() - 1] << endl;

  double w = (4 * nZ + 2 * lZ - 1) / 5.;
  double V0 = fBetaType * 3 * alpha * Zd / 2. / R;
  // double V0 = -alpha*sqrt(2.*1.825E-4/M_PI)*26.;
  double e = (sqr(W0 - W) + sqr(W + V0) - 1) / 6.;

  double Ap = 1.;
  double sum = 0.;
  for (int j = 0; j < occNumbersZ.size(); j += 4) {
    if (occNumbersZ[j + 1] == 0) {
      sum += occNumbersZ[j + 3];
    }
  }
  Ap = (2. * (Zd - fBetaType) / sum - 2.) / 3.;

  // cout << "Ap: " << Ap << endl;

  return 1 - 8. / 5. * w * e * R * R / (5. * Ap + 2);

  double weakR2;
  if (An > 90) {
    nu = utilities::CalcChargeIndepNu(rms, Zd - fBetaType,
                                      An - (Zd - fBetaType));
    weakR2 = utilities::WeakIntegratedRMS(nu, nN, lN, nu, nZ, lZ);
  } else {
    nu = utilities::CalcNu(rms, Zd - fBetaType);
    if (fBetaType == BETA_MINUS) {
      weakR2 = pow(utilities::GetRMSHO(nN, lN, nu), 2);
    } else {
      weakR2 = pow(utilities::GetRMSHO(nZ, lZ, nu), 2);
    }
  }

  /*cout << "Nu: " << nu << endl;
  cout << "nN: " << nN << " lN: " << lN << endl;
  cout << "nZ: " << nZ << " lZ: " << lZ << endl;*/

  double delta = (weakR2 - 3. / 5. * R * R) / 6.;

  // return (1-e*weakR2)/(1-e*3./5.*R*R);
  // return 1-e*(weakR2-R*R*3/5.);

  // cout << "rmsW: " << weakR2 << " rmsCh " << 3./5.*R*R << " Delta: " << delta
  // << endl;

  AC0 = alpha * Zd / R * (fBetaType * 4. / 3. * W0 - 9 / 2. * alpha * Zd / R +
                          4 * alpha * Zd / R * W0 * W0 * delta);
  AC1 = 40 / 9. * W0 +
        alpha * Zd / R *
            (fBetaType * 16. / 3. * W0 * W0 * delta - fBetaType * 22. / 3. -
             22 * alpha * Zd / R * W0 * delta);
  AC2 = -40 / 9. - fBetaType * 24 * alpha * Zd / R * W0 * delta +
        27 * pow(alpha * Zd / R, 2) * delta + 44. / 9. * W0 * W0 * delta;
  ACm1 = -4 / 9. * W0;

  /*cout << "W: " << W << " W0: " << W0 << endl;
  cout << "AC0: " << AC0 << " AC1: " << AC1 << " AC2: " << AC2 << " ACm1: " <<
  ACm1 << endl;*/

  VC0 = alpha * Zd / R * (-fBetaType * 4 * W0 - 9. / 2. * alpha * Zd / R +
                          4 * alpha * Zd / R * W0 * W0 * delta);
  VC1 = 8. / 3. * W0 +
        alpha * Zd / R * (-fBetaType * 16. / 3. * W0 * W0 * delta -
                          fBetaType * 2 - 2. * alpha * Zd / R * W0 * delta);
  VC2 = -8 / 3. + fBetaType * 8. / 3. * alpha * Zd / R * W0 * delta +
        7. * pow(alpha * Zd / R, 2) * delta;
  VCm1 = 4. / 3. * W0;

  double r2 = 0.75 * R + weakR2 / 2. / R;
  double r2r = 0.7 * weakR2 / R;

  /*VC0 = -alpha*Zd*2./9.*W0*(r2r/2.-0.25*R);
  VC1 = alpha*Zd*(4./3.*r2+4./9.*r2r-14./9.*R);
  VCm1 = alpha*Zd*(2*r2-r2r-0.5*R);
  VC2 = 0.;*/

  if (fDecayType == FERMI) {
    return 1 + delta * (VC0 + VC1 * W + VC2 * W * W + VCm1 / W);
  } else if (fDecayType == GAMOW_TELLER) {
    return 1 + delta * (AC0 + AC1 * W + AC2 * W * W + ACm1 / W);
  } else if (MixingRatio > 0) {
    return 1 +
           delta * (1. / (1 + pow(MixingRatio, 2)) *
                        (VC0 + VC1 * W + VC2 * W * W + VCm1 / W) +
                    1 / (1 + 1. / pow(MixingRatio, 2)) *
                        (AC0 + AC1 * W + AC2 * W * W + ACm1 / W));
  }
  cout << "Mixing ratio is badly defined. Returning 1." << endl;

  return 1.;
}

double SpectralFunctions::RelativisticCorrection(double W) {
  double Wb = W + fBetaType * 3 * alpha * Zd / (2. * R);
  double pb = sqrt(Wb * Wb - 1.);
  double H2 = -pow(pb * R, 2.) / 6.;
  double D1 = Wb * R / 3.;
  double D3 = -Wb * R * pow(pb * R, 2.) / 30 - fBetaType * alpha * Zd / 10.;
  double d1 = R / 3.;
  double d3 = -R * pow(pb * R, 2.) / 30.;
  double N1 = (W0 - W) * R / 3.;
  double N2 = -pow((W0 - W) * R, 2.) / 6.;
  double N3 = pow((W0 - W) * R, 3.) / 30;

  double Vf2, Vf3;
  double Af2, Af3;

  Vf2 = -2. * (D1 + N1) + 2 * gamma / W * d1;
  Vf3 = -2. * (D3 + N1 * H2 - N2 * D1 - N3) + 2 * gamma / W * (d3 - N2 * d1);

  Af2 = 2. / sqrt(3) * (D1 + N1) - 2. / sqrt(3) * gamma / W * d1;
  Af3 = 2. * sqrt(2. / 3.) * (D1 - N1) - 2. * sqrt(2. / 3.) * gamma / W * d1;

  double mismatch = W0 - fBetaType * 2.5 + fBetaType * 6. / 5. * alpha * Zd / R;

  if (fDecayType == FERMI) {
    return 1. - 3. / 10 * R * mismatch * Vf2 - 3. / 28. * R * mismatch * Vf3;
  }

  vector<int> occNumbersZ = utilities::GetOccupationNumbers(Zd - fBetaType);
  vector<int> occNumbersN =
      utilities::GetOccupationNumbers(An - Zd + fBetaType);

  int li, lf;
  return 1;
}

double SpectralFunctions::DeformationCorrection(double W) {
  double bOverA = utilities::CalcBoverA(beta2);
  double a, b;

  a = R * sqrt(3. / (bOverA * bOverA + 2.));
  b = bOverA * a;

  // cout << "b/a, a, b" << bOverA << " " << a << " " << b <<endl;

  double V0s = 3. * alpha * Zd / 2. / R;
  double V0d = 0.;

  if (beta2 > 0) {
    V0d = 3. * alpha * Zd / 2. / a / sqrt(bOverA * bOverA - 1.) *
          log(sqrt(bOverA * bOverA - 1) + bOverA);
  } else if (beta2 < 0) {
    V0d = 3. * alpha * Zd / 2. / a / sqrt(1. - bOverA * bOverA) *
          asin(sqrt(1 - bOverA * bOverA));
  }
  double etaS = 1. / 6. * (pow(W0 - W, 2.) + pow(W + fBetaType * V0s, 2.) - 1.);
  double etaD = 1. / 6. * (pow(W0 - W, 2.) + pow(W + fBetaType * V0d, 2.) - 1.);
  double DC0 = (1 - etaD * 3. / 5. * R * R) / (1 - etaS * 3. / 5. * R * R);

  // cout << "DC0 " << DC0 << endl;

  int steps = 1E2;
  double grid[steps];
  double integrand[steps];

  double prefact = 4. / 3. * M_PI * pow(R, 2. * (1 - gamma));

  for (int i = 0; i < steps; i++) {
    grid[i] = a + (i + 1.) / (double)steps * (b - a);
    integrand[i] = pow(grid[i], 3.) * Zd *
                   abs(utilities::GetDerivDeformedChargeDist(grid[i], a, b)) *
                   L0Correction(W, grid[i]) * pow(grid[i], 2. * (gamma - 1));
    cout << grid[i] << "\t" << integrand[i] << endl;
  }

  double newL0 = prefact * utilities::Trapezoid(grid, integrand, steps);

  // cout << "newL0 " << newL0 << endl;
  double DFS = newL0 / L0Correction(W, R);

  DC0 = 1.;

  return DC0 * DFS;
}

double SpectralFunctions::L0Correction(double W, double r) {
  double sum = 0;
  double common = 0;
  double specific = 0;
  for (int i = 1; i < 7; i++) {
    if (fBetaType == BETA_PLUS)
      sum += aPos[i] * pow(W * r, i - 1);
    else
      sum += aNeg[i] * pow(W * r, i - 1);
  }
  common = (1 + 13. / 60. * pow(alpha * Zd, 2) -
            fBetaType * W * r * alpha * Zd * (41. - 26. * gamma) / 15. /
                (2. * gamma - 1) -
            fBetaType * alpha * Zd * r * gamma * (17. - 2. * gamma) / 30. / W /
                (2. * gamma - 1) +
            sum);
  if (fBetaType == BETA_PLUS)
    specific = aPos[0] * r / W + 0.22 * (r - 0.0164) * pow(alpha * Zd, 4.5);
  else
    specific = aNeg[0] * r / W + 0.41 * (r - 0.0164) * pow(alpha * Zd, 4.5);
  return (common + specific) * 2 / (1 + gamma);
}

double SpectralFunctions::UCorrection(double W) {
  double a0 = -5.6E-5 - fBetaType * 4.94E-5 * Zd + 6.23E-8 * pow(Zd, 2);
  double a1 = 5.17E-6 + fBetaType * 2.517E-6 * Zd + 2.00E-8 * pow(Zd, 2);
  double a2 = -9.17e-8 + fBetaType * 5.53E-9 * Zd + 1.25E-10 * pow(Zd, 2);

  double p = sqrt(W * W - 1);

  return 1 + a0 + a1 * p + a2 * p * p;
}

double SpectralFunctions::QCorrection(double W) {
  double a = 0;

  if (fDecayType == FERMI)
    a = 1.;
  else if (fDecayType == GAMOW_TELLER)
    a = -1. / 3.;
  else if (MixingRatio > 0)
    a = (1 - pow(MixingRatio, 2.) / 3) / (1 + pow(MixingRatio, 2));

  double M = An * (protonMasskeV + neutronMasskeV) / 2. / electronMasskeV;

  double p = sqrt(W * W - 1);

  return 1 - fBetaType * M_PI * alpha * Zd / M / p * (1 + a * (W0 - W) / 3 / M);
}

double SpectralFunctions::RadiativeCorrection(double W) {
  // 1st order, based on the 5th Wilkinson article
  double beta = sqrt(1.0 - 1.0 / W / W);

  double g = 3 * log(protonMasskeV / electronMasskeV) - 0.75 +
             4 * (atanh(beta) / beta - 1) *
                 ((W0 - W) / 3 / W - 1.5 + log(2 * (W0 - W)));
  g += 4.0 / beta * Spence(2 * beta / (1 + beta)) +
       atanh(beta) / beta * (2 * (1 + beta * beta) +
                             (W0 - W) * (W0 - W) / 6 / W / W - 4 * atanh(beta));

  double O1corr = alpha / 2 / M_PI *
                  (g - 3 * log(protonMasskeV / electronMasskeV / 2 / W0));

  double L = 1.026725 * pow(1 - 2 * alpha / 3 / M_PI * log(2 * W0), 9 / 4.);

  // 2nd order
  double d1f, d2, d3, d14;
  double Lambda = sqrt(10) / R;
  double LambdaOverM =
      Lambda / protonMasskeV * electronMasskeV;  // this is dimensionless
  double gamma_E = 0.5772;

  d14 =
      log(protonMasskeV / electronMasskeV) - 5. / 3. * log(2 * W) + 43. / 18.;

  d1f = log(LambdaOverM) - gamma_E + 4. / 3. - log(sqrt(10.0)) -
        3.0 / M_PI / sqrt(10.0) * (0.5 + gamma_E + log(sqrt(10) / LambdaOverM));

  d2 = 3.0 / 2.0 / M_PI / sqrt(10.0) * LambdaOverM *
       (1 - M_PI / 2 / sqrt(10) * LambdaOverM);

  d3 = 3.0 * g_A * g_M / M_PI / sqrt(10.0) * LambdaOverM *
       (gamma_E - 1 + log(sqrt(10) / LambdaOverM) +
        M_PI / 4 / sqrt(10) * LambdaOverM);

  double O2corr = fBetaType * alpha * alpha * Zd * (d14 + d1f + d2 + d3);

  // 3rd order
  double a = 0.5697;
  double b = 4 / 3 / M_PI * (11 / 4. - gamma_E - M_PI * M_PI / 6);
  double f = log(2 * W0) - 5. / 6.;
  double g2 =
      0.5 * (pow(log(R), 2.) - pow(log(2 * W), 2.)) + 5. / 3. * log(2 * R * W);

  double O3corr =
      pow(alpha, 3) * pow(Zd, 2) *
      (a * log(Lambda / W) + b * f + 4 * M_PI / 3. * g2 - 0.649 * log(2 * W0));

  return (1 + O1corr) * (L + O2corr + O3corr);
}

double SpectralFunctions::NeutrinoRadiativeCorrection(double Wv) {
  double h = 0.;
  double pv = sqrt(Wv * Wv - 1);
  double beta = pv / Wv;
  h += 3 * log(protonMasskeV / electronMasskeV) + 23 / 4. +
       8 / beta * Spence(2 * beta / (1 + beta)) +
       8 * (atan(beta) / beta - 1) * log(2 * Wv * beta) +
       4 * atan(beta) / beta * ((7 + 3 * beta * beta) / 8 - 2 * atan(beta));
  return 1 + alpha / 2 / M_PI * h;
}

double SpectralFunctions::Spence(double x) { return -gsl_sf_dilog(x); }

double SpectralFunctions::RecoilCorrection(double W) {
  double Vr0, Vr1, Vr2, Vr3;
  double Ar0, Ar1, Ar2, Ar3;
  double MassOfNucleus = An * (protonMasskeV + neutronMasskeV) / 2. /
                         electronMasskeV;  // in units of electron mass
  double MassOfNucleus2 = sqr(MassOfNucleus);

  Ar0 = -2. * W0 / 3. / MassOfNucleus - W0 * W0 / 6. / MassOfNucleus2 -
        77. / 18. / MassOfNucleus2;
  Ar1 = -2. / 3. / MassOfNucleus + 7. * W0 / 9. / MassOfNucleus2;
  Ar2 = 10. / 3. / MassOfNucleus - 28. * W0 / 9. / MassOfNucleus2;
  Ar3 = 88. / 9. / MassOfNucleus2;

  Vr0 = W0 * W0 / 2. / MassOfNucleus2 - 11. / 6. / MassOfNucleus2;
  Vr1 = W0 / 3. / MassOfNucleus2;
  Vr2 = 2. / MassOfNucleus - 4. * W0 / 3. / MassOfNucleus2;
  Vr3 = 16. / 3. / MassOfNucleus2;

  if (fDecayType == FERMI) {
    return 1 + Vr0 + Vr1 / W + Vr2 * W + Vr3 * W * W;
  } else if (fDecayType == GAMOW_TELLER) {
    return 1 + Ar0 + Ar1 / W + Ar2 * W + Ar3 * W * W;
  } else if (MixingRatio > 0) {
    return 1 +
           1. / (1 + pow(MixingRatio, 2)) *
               (Vr0 + Vr1 / W + Vr2 * W + Vr3 * W * W) +
           1. / (1 + 1. / pow(MixingRatio, 2)) *
               (Ar0 + Ar1 / W + Ar2 * W + Ar3 * W * W);
  }
  cout << "Mixing ratio badly defined. Returning 1." << endl;
  return 1;
}

double SpectralFunctions::AtomicScreeningCorrection(double W) {
  vector<double> Aby, Bby;

  utilities::PotParam(Zd - 1 * fBetaType, Aby, Bby);

  double l = 2 * (Aby[0] * Bby[0] + Aby[1] * Bby[1] + Aby[2] * Bby[2]);

  double p = sqrt(W * W - 1);

  double Wt = W - fBetaType * 0.5 * alpha * (Zd - fBetaType) * l;

  complex<double> pt;

  pt = 0.5 * p +
       0.5 * sqrt(complex<double>(p * p - fBetaType * 2 * alpha * Zd * Wt * l));

  double y = fBetaType * alpha * Zd * W / p;
  complex<double> yt = fBetaType * alpha * Zd * Wt / pt;

  /*cout << "V: " << W-Wt << endl;
  cout << W << " " << yt.real() << " " << yt.imag() << endl;
  cout << W << " " << pt.real() << " " << pt.imag() << endl;*/

  gsl_sf_result magn;
  gsl_sf_result phase;
  gsl_sf_lngamma_complex_e(gamma, y, &magn, &phase);

  gsl_sf_result magnT;
  gsl_sf_result phaseT;
  gsl_sf_lngamma_complex_e(gamma - yt.imag(), yt.real(), &magnT, &phaseT);
  // now we have what we want in magn.val

  double first = Wt / W;
  double second = exp(2 * (magnT.val - magn.val));

  gsl_sf_lngamma_complex_e(gamma - 2 * pt.imag() / l, 2 * pt.real() / l, &magnT,
                           &phaseT);
  gsl_sf_lngamma_complex_e(1, 2 * p / l, &magn, &phase);
  double third = exp(2 * (magnT.val - magn.val));
  double fourth = exp(-M_PI * y);
  double fifth = pow(2 * p / l, 2 * (1 - gamma));

  return first * second * third * fourth * fifth;
}

double SpectralFunctions::AtomicExchangeCorrection(double W) {
  double E = W - 1;

  return 1 + exPars[0] / E + exPars[1] / E / E +
         exPars[2] * exp(-exPars[3] * E) +
         exPars[4] * sin(pow(W - exPars[6], exPars[5]) + exPars[7]) /
             pow(W, exPars[8]);
}

double SpectralFunctions::AtomicMismatchCorrection(double W) {
  double dBdZ2 = (44.200 * pow(Zd - fBetaType, 0.41) +
                  2.3196E-7 * pow(Zd - fBetaType, 4.45)) /
                 electronMasskeV / 1000.;

  double K = -0.872 + 1.270 * pow(Zd, 0.097) + 9.062E-11 * pow(Zd, 4.5);
  double vp = sqrt(1 - 1 / W / W);
  double l = 1.83E-3 * K * Zd / vp;
  double M = An * (protonMasskeV + neutronMasskeV) / 2. / electronMasskeV;
  // assume as an average the recoil velocity at half-momentum transfer ~
  // sqrt(W0^2-1)/2
  double vR = sqrt(1 - M * M / (M * M + (W0 * W0 - 1) / 4.));

  double psi2 = 1 + 2 * alpha / vp * (atan(1 / l) - l / 2 / (1 + l * l));

  double C0 = -alpha * alpha * Zd * alpha / vp * l / (1 + l * l) / psi2;

  double C1 = 2 * alpha * alpha * Zd * vR / vp *
              ((0.5 + l * l) / (1 + l * l) - l * atan(1 / l)) / psi2;

  return 1 - 2 / (W0 - W) * (0.5 * dBdZ2 + 2 * (C0 + C1));
}

double SpectralFunctions::CalculateWeakMagnetism() {
  utilities::NuclearState nsf = nilsson::CalculateDeformedState(Zd, An, beta2, beta4);
  utilities::NuclearState nsi =
      nilsson::CalculateDeformedState(Zd + fBetaType, An, beta2, beta4);

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

vector<double> SpectralFunctions::GetSpectrumShape(double Energy) {
  double W = Energy / electronMasskeV + 1;

  double result = 1;
  double neutrinoResult = 1;

  double Wv = W0 - W + 1;

  if (fPhaseSpace) {
    result *= sqrt(W * W - 1) * W * sqr(W0 - W);
    neutrinoResult *= sqrt(Wv * Wv - 1) * Wv * sqr(W0 - Wv);
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

  if (fAtomicExchangeCorrection) {
    if (fBetaType == BETA_MINUS) {
      result *= AtomicExchangeCorrection(W);
      neutrinoResult *= AtomicExchangeCorrection(Wv);
    } else {
      // cout << "The exchange correction does not exist for beta+ decay." <<
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
  cout << "Warning: Deformation corrections and relativistic corrections for "
          "GT transitions have not yet been implemented." << endl;
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
