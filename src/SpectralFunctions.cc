#include "BSG/SpectralFunctions.h"
#include "BSG/Screening.h"

#include "PDS/Units/GlobalSystemOfUnits.h"
#include "PDS/Units/GlobalPhysicalConstants.h"
#include "NHL/Units/BetaDecayPhysicalConstants.h"
#include "NHL/ShellModel.h"
#include "NHL/ChargeDistributions.h"

#include <complex>

// GNU Scientific Library stuff
// http://www.gnu.org/software/gsl/
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_const_mksa.h"
#include "gsl/gsl_const_num.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_result.h"
#include "gsl/gsl_sf_dilog.h"
#include "gsl/gsl_integration.h"

double BSG::SpectralFunctions::PhaseSpace(double W, double W0) {
  double result = std::sqrt(W * W - 1.) * W * std::pow(W0 - W, 2.);
  return result;
}

double BSG::SpectralFunctions::FermiFunction(double W, int Z, double R,
                                        int betaType) {
  if (Z == 0) {
    return 1.;
  }
  double gamma = std::sqrt(1. - std::pow(fine_structure_const * Z, 2.));
  double p = std::sqrt(W * W - 1.);
  double y = betaType * fine_structure_const * Z * W / p;
  //IMPORTANT
  //Note that we use here 'traditional' definition of the Fermi function
  //using 4 instead of 2 * (gamma + 1) as a prefactor
  //See Reviews of Modern Physics 90 (2018) 015008
  double first = 4.;
  // the second term will be incorporated in the fifth
  // double second = 1/std::pow(gsl_sf_gamma(2*gamma+1),2);
  double second = std::pow(2. * p * R, 2. * (gamma - 1.));
  double third = std::exp(pi * y);

  // the fourth is a bit tricky
  // we use the complex gamma function from GSL
  gsl_sf_result magn;
  gsl_sf_result phase;
  gsl_sf_lngamma_complex_e(gamma, y, &magn, &phase);
  // now we have what we wAt in magn.val

  // but we incorporate the other term here as well
  double fourth = std::exp(2. * (magn.val - gsl_sf_lngamma(2. * gamma + 1.)));

  double result = first * second * third * fourth;
  return result;
}

double BSG::SpectralFunctions::L0Correction(double W, int Z, double r, int betaType,
                                       std::array<double, 7> aPos, std::array<double, 7> aNeg) {
  double gamma = std::sqrt(1. - std::pow(fine_structure_const * Z, 2.));
  double sum = 0;
  double common = 0;
  double specific = 0;
  for (unsigned int i = 1; i < aPos.size(); i++) {
    if (betaType == NHL::BETA_PLUS)
      sum += aPos[i] * std::pow(W * r, i - 1);
    else
      sum += aNeg[i] * std::pow(W * r, i - 1);
  }
  common = 1. + 13. / 60. * std::pow(fine_structure_const * Z, 2) -
           betaType * W * r * fine_structure_const * Z * (41. - 26. * gamma) / 15. /
               (2. * gamma - 1) -
           betaType * fine_structure_const * Z * r * gamma * (17. - 2. * gamma) / 30. / W /
               (2. * gamma - 1) +
           sum;
  if (betaType == NHL::BETA_PLUS)
    specific = aPos[0] * r / W + 0.22 * (r - 0.0164) * std::pow(fine_structure_const * Z, 4.5);
  else
    specific = aNeg[0] * r / W + 0.41 * (r - 0.0164) * std::pow(fine_structure_const * Z, 4.5);
  return common + specific;
}

double BSG::SpectralFunctions::UCorrection(double W, int Z, double R, int betaType,
                                      NHL::NuclearShapes ESShape,
                                      std::vector<double>& v,
                                      std::vector<double>& vp) {
  double result = 1.;
  if (ESShape == NHL::NuclearShapes::FERMI) {
    double a0 = -5.6E-5 - betaType * 4.94E-5 * Z + 6.23E-8 * std::pow(Z, 2);
    double a1 = 5.17E-6 + betaType * 2.517E-6 * Z + 2.00E-8 * std::pow(Z, 2);
    double a2 = -9.17e-8 + betaType * 5.53E-9 * Z + 1.25E-10 * std::pow(Z, 2);

    double p = std::sqrt(W * W - 1);

    result = 1. + a0 + a1 * p + a2 * p * p;
  }
  //TODO check
  else {
    result *= UCorrection(W, Z, R, betaType, v, vp);
  }

  return result;
}

double BSG::SpectralFunctions::UCorrection(double W, int Z, double R, int betaType,
                                      std::vector<double>& v,
                                      std::vector<double>& vp) {
  double delta1 = 4. / 3. * (vp[0] - v[0]) + 17. / 30. * (vp[1] - v[1]) +
                  25. / 63. * (vp[2] - v[2]);
  double delta2 = 2. / 3. * (vp[0] - v[0]) + 7. / 12. * (vp[1] - v[1]) +
                  11. / 6. * (vp[2] - v[2]);
  double delta3 = 1. / 3. * (vp[0] * vp[0] - v[0] * v[0]) +
                  1. / 15. * (vp[1] * vp[1] - v[1] * v[1]) +
                  1. / 35. * (vp[2] * vp[2] - v[2] * v[2]) +
                  1. / 6. * (vp[1] * vp[0] - v[1] * v[0]) +
                  1. / 9. * (vp[2] * vp[0] - v[2] * v[0]) +
                  1. / 20. * (vp[2] * vp[1] - v[2] * v[1]) +
                  1. / 5. * (vp[1] - v[1]) + 1. / 7. * (vp[2] - v[2]);
  double delta4 = 4. / 3. * (vp[0] - v[0]) + 4. / 5. * (vp[1] - v[1]) +
                  4. / 7. * (vp[2] - v[2]);

  double gamma = std::sqrt(1. - (fine_structure_const * Z) * (fine_structure_const * Z));

  double result = 1. + betaType * fine_structure_const * Z * W * R * delta1 +
                  betaType * gamma / W * fine_structure_const * Z * R * delta2 +
                  (fine_structure_const * Z) * (fine_structure_const * Z) * delta3 -
                  (W * R) * (W * R) * delta4;

  return result;
}

double BSG::SpectralFunctions::QCorrection(double W, double W0, int Z, int A,
                                      int betaType, NHL::BetaDecayType decayType,
                                      double mixingRatio) {
  double a = 0;

  if (decayType == NHL::BetaDecayType::FERMI)
    a = 1.;
  else if (decayType == NHL::BetaDecayType::GAMOWTELLER)
    a = -1. / 3.;
  else if (mixingRatio > 0.)
    a = (1. - std::pow(mixingRatio, 2.) / 3.) / (1. + std::pow(mixingRatio, 2));

  double M = A * NHL::nucleon_mass_bu;

  double p = std::sqrt(W * W - 1.);

  return 1. - betaType * pi * fine_structure_const * Z / M / p * (1. + a * (W0 - W) / 3. / M);
}

//TODO Check for beta+
double BSG::SpectralFunctions::RadiativeCorrection(double W, double W0, int Z,
                                              double R, int betaType, double gA,
                                              double gM) {
  // 1st order, based on the 5th Wilkinson article
  double beta = std::sqrt(1.0 - 1.0 / W / W);

  double g = 3. * std::log(proton_mass_c2 / electron_mass_c2) - 0.75 +
             4. * (std::atanh(beta) / beta - 1.) *
                 ((W0 - W) / 3. / W - 1.5 + std::log(2 * (W0 - W)));
  g += 4.0 / beta * Spence(2. * beta / (1. + beta)) +
       std::atanh(beta) / beta *
           (2. * (1. + beta * beta) + (W0 - W) * (W0 - W) / 6. / W / W -
            4. * std::atanh(beta));

  double O1corr =
      fine_structure_const / 2. / pi *
      (g - 3. * std::log(proton_mass_c2 / electron_mass_c2 / 2. / W0));

  double L =
      1.026725 * std::pow(1. - 2. * fine_structure_const / 3. / pi * std::log(2. * W0), 9. / 4.);

  // 2nd order
  double d1f, d2, d3, d14;
  double lambda = std::sqrt(10) / R;
  double lambdaOverM =
      lambda / NHL::nucleon_mass_bu;  // this is dimensionless

  d14 = std::log(proton_mass_c2 / electron_mass_c2) -
        5. / 3. * std::log(2 * W) + 43. / 18.;

  d1f = std::log(lambdaOverM) - euler_mascheroni_constant + 4. / 3. -
        std::log(std::sqrt(10.0)) -
        3.0 / pi / std::sqrt(10.0) * lambdaOverM * (0.5 + euler_mascheroni_constant +
                                        std::log(std::sqrt(10) / lambdaOverM));

  d2 = 3.0 / 2.0 / pi / std::sqrt(10.0) * lambdaOverM *
       (1. - pi / 2. / std::sqrt(10) * lambdaOverM);

  d3 = 3.0 * gA * gM / pi / std::sqrt(10.0) * lambdaOverM *
       (euler_mascheroni_constant - 1. + std::log(std::sqrt(10) / lambdaOverM) +
        pi / 4 / std::sqrt(10) * lambdaOverM);

  double O2corr = fine_structure_const * fine_structure_const * Z * (d14 + d1f + d2 + d3);

  // 3rd order
  double a = 0.5697;
  double b =
      4. / 3. / pi * (11. / 4. - euler_mascheroni_constant - pi * pi / 6);
  double f = std::log(2 * W) - 5. / 6.;
  double g2 = 0.5 * (std::pow(std::log(R), 2.) - std::pow(std::log(2 * W), 2.)) +
              5. / 3. * std::log(2 * R * W);

  double O3corr = std::pow(fine_structure_const, 3) * std::pow(Z, 2) *
                  (a * std::log(lambda / W) + b * f + 4. / pi / 3. * g2 -
                   0.649 * std::log(2 * W0));

  return (1 + O1corr) * (L + O2corr + O3corr);
}

double BSG::SpectralFunctions::RecoilCorrection(double W, double W0, int A,
                                           NHL::BetaDecayType decayType, double mixingRatio) {
  double Vr0, Vr1, Vr2, Vr3;
  double Ar0, Ar1, Ar2, Ar3;
  double M = A * NHL::nucleon_mass_bu;  // in units of electron mass
  double M2 = std::pow(M, 2.);

  Ar0 = -2. * W0 / 3. / M - W0 * W0 / 6. / M2 - 77. / 18. / M2;
  Ar1 = -2. / 3. / M + 7. * W0 / 9. / M2;
  Ar2 = 10. / 3. / M - 28. * W0 / 9. / M2;
  Ar3 = 88. / 9. / M2;

  Vr0 = W0 * W0 / 2. / M2 - 11. / 6. / M2;
  Vr1 = W0 / 3. / M2;
  Vr2 = 2. / M - 4. * W0 / 3. / M2;
  Vr3 = 16. / 3. / M2;

  if (decayType == NHL::BetaDecayType::FERMI) {
    return 1 + Vr0 + Vr1 / W + Vr2 * W + Vr3 * W * W;
  } else if (decayType == NHL::BetaDecayType::GAMOWTELLER) {
    return 1 + Ar0 + Ar1 / W + Ar2 * W + Ar3 * W * W;
  } else if (mixingRatio > 0) {
    return 1 +
           1. / (1 + std::pow(mixingRatio, 2)) *
               (Vr0 + Vr1 / W + Vr2 * W + Vr3 * W * W) +
           1. / (1 + 1. / std::pow(mixingRatio, 2)) *
               (Ar0 + Ar1 / W + Ar2 * W + Ar3 * W * W);
  }
  return 1;
}

std::tuple<double, double> BSG::SpectralFunctions::CCorrectionComponents(
    double W, double W0, int Z, double R, int betaType, NHL::BetaDecayType decayType,
    double gA, double gP, double bAc, double dAc, double lambda,
    NHL::NuclearShapes NSShape, double modGaussFit) {

    double AC0, AC1, ACm1, AC2;
    double VC0, VC1, VCm1, VC2;

    //Uniformly charged sphere results
    double F1111 = 27./35.;
    double F1221 = 57./70.;
    double F1222 = 233./210.;
    double F1211 = -3./70.;

    if (NSShape == NHL::NuclearShapes::MODGAUSS) {
      F1111 = 0.757 + 0.0069 * (1 - std::exp(-modGaussFit / 1.008));
      F1221 = 0.844 - 0.0182 * (1 - std::exp(-modGaussFit / 1.974));
      F1222 = 1.219 - 0.0640 * (1 - std::exp(-modGaussFit / 1.550));
    }

    VC0 = -std::pow(W0 * R, 2.) / 5. -
          betaType * 2. / 9. * fine_structure_const * Z * W0 * R * F1111 -
          std::pow(fine_structure_const * Z, 2.) / 3. * F1222;

    VC1 = 4. / 15. * W0 * R * R -
          betaType * 2. / 3. * fine_structure_const * Z * R * (F1221 - F1111 / 3.);

    VCm1 = 2. / 15. * W0 * R * R - betaType * fine_structure_const * Z * R / 3. * F1211;

    VC2 = -4. / 15. * R * R;

    AC0 = -1. / 3. * std::pow(fine_structure_const * Z, 2.) * F1222 -
          1. / 5. * (W0 * W0 - 1.) * R * R +
          betaType * 2. / 27. * fine_structure_const * Z * W0 * R * F1111 + 11. / 45. * R * R;

    AC1 = 4. / 9. * W0 * R * R -
          betaType * 2. / 3. * fine_structure_const * Z * R * (1. / 9. * F1111 + F1221);

    ACm1 = -2. / 45. * W0 * R * R + betaType * fine_structure_const * Z * R / 3. * F1211;

    AC2 = -4. / 9. * R * R;

    double cShape = 0.;

    if (decayType == NHL::BetaDecayType::FERMI) {
      cShape = 1. + VC0 + VC1 * W + VCm1 / W + VC2 * W * W;
    } else if (decayType == NHL::BetaDecayType::GAMOWTELLER) {
      cShape = 1. + AC0 + AC1 * W + ACm1 / W + AC2 * W * W;
    }

    double cNS = 0;
    if (decayType == NHL::BetaDecayType::GAMOWTELLER) {
      double Mn = NHL::nucleon_mass_bu;

      //double Lambda = std::sqrt(2.)/3.*10.*ratioM121;

      double phi = gP / gA / std::pow(2.* Mn * R, 2.);

      double NSC0 = -1. / 45. * R * R * lambda +
                    1. / 3. * W0 / Mn * (-betaType * 2. * bAc + dAc) +
                    betaType * 2. / 5. * fine_structure_const * Z / Mn / R *
                        (betaType * 2. * bAc + dAc) -
                    betaType * 2. / 35. * fine_structure_const * Z * W0 * R * lambda;

      double NSC1 = betaType * 4. / 3. / Mn * bAc  -
                    2. / 45. * W0 * R * R * lambda +
                    betaType * fine_structure_const * Z * R * 2. / 35. * lambda;

      double NSCm1 = -1. / 3. / Mn * (betaType * 2. * bAc + dAc) +
                     2. / 45. * W0 * R * R * lambda;

      double NSC2 = 2. / 45. * R * R * lambda;

      double gamma = std::sqrt(1. - std::pow(fine_structure_const * Z, 2.));

      double P0 = betaType * 2. / 25. * fine_structure_const * Z * R * W0 + 51. / 250. * std::pow(fine_structure_const*Z, 2.);
      double P1 = betaType * 2. / 25. * fine_structure_const * Z * R;
      double Pm1 = -2. /3. * gamma * W0 * R * R + betaType * 26. / 25. * fine_structure_const * Z * R * gamma;

      cNS = NSC0 + NSC1 * W + NSCm1 / W + NSC2 * W * W;

      cNS += phi * (P0 + P1 * W + Pm1 / W);
    }

    return std::make_tuple(cShape, cNS);
}

double BSG::SpectralFunctions::CCorrection(double W, double W0, int Z, int A,
                                      double R, int betaType,
                                      NHL::BetaDecayType decayType, double gA, double gP,
                                      double bAc, double dAc, double lambda, bool addCI,
                                      NHL::NuclearShapes NSShape, double hoFit) {
  double cShape, cNS;
  std::tie(cShape, cNS) =
      CCorrectionComponents(W, W0, Z, R, betaType, decayType, gA, gP,
                            bAc, dAc, lambda, NSShape, hoFit);
  double result = 0.;
  if (addCI) {
    result = cShape * CICorrection(W, W0, Z, A, R, betaType) + cNS;
  } else {
    result = cShape + cNS;
  }
  return result;
}

double BSG::SpectralFunctions::CCorrection(
    double W, double W0, int Z, double R, int betaType,
    NHL::BetaDecayType decayType, double gA, double gP, double bAc, double dAc,
    double lambda, bool addCI, NHL::NuclearShapes NSShape, double hoFit,
    NHL::SingleParticleState& spsi,
    NHL::SingleParticleState& spsf) {

  double cShape, cNS;
  std::tie(cShape, cNS) =
      CCorrectionComponents(W, W0, Z, R, betaType, decayType, gA, gP,
                            bAc, dAc, lambda, NSShape, hoFit);
  double result = 0.;
  if (addCI) {
    result = cShape * CICorrection(W, W0, Z, R, betaType, spsi, spsf) + cNS;
  } else {
    result = cShape + cNS;
  }
  return result;
}

double BSG::SpectralFunctions::CICorrection(double W, double W0, int Z, int A,
                                       double R, int betaType) {
  /*double rms = std::sqrt(3. / 5.) * R;
  double nu = 0.;*/

  //TODO check

  int nN, lN, nZ, lZ;
  std::vector<int> occNumbersN =
      NHL::GetOccupationNumbers(A - (Z - betaType));
  nN = occNumbersN[occNumbersN.size() - 1 - 3];
  lN = occNumbersN[occNumbersN.size() - 1 - 2];
  std::vector<int> occNumbersZ = NHL::GetOccupationNumbers(Z - betaType);
  nZ = occNumbersZ[occNumbersZ.size() - 1 - 3];
  lZ = occNumbersZ[occNumbersZ.size() - 1 - 2];

  double w = (4 * nZ + 2 * lZ - 1) / 5.;
  double V0 = betaType * 3 * fine_structure_const * Z / 2. / R;
  double e = (std::pow(W0 - W, 2.) + std::pow(W + V0, 2.) - 1) / 6.;

  double Ap = 1.;
  double sum = 0.;
  for (unsigned int j = 0; j < occNumbersZ.size(); j += 4) {
    if (occNumbersZ[j + 1] == 0) {
      sum += occNumbersZ[j + 3];
    }
  }
  Ap = (2. * (Z - betaType) / sum - 2.) / 3.;

  return 1 - 8. / 5. * w * e * R * R / (5. * Ap + 2);
}

//TODO
double BSG::SpectralFunctions::CICorrection(
    double W, double W0, int Z, double R, int betaType,
    NHL::SingleParticleState& spsi,
    NHL::SingleParticleState& spsf) {
/*  double V0 = betaType * 3. * Z * fine_structure_const / 2. / R;
  double epsilon = 1. / 6. * (std::pow(W0 - W2., ) + std::pow(W + V0, 2.) - 1.);

  double nu = NHL::CalcNu(R * std::sqrt(3. / 5.), Z);

  double result = 0.;

  double C = 0.;
  for (int i = 0; i < spsi.componentsHO.size(); i++) {
    for (int j = 0; j < spsf.componentsHO.size(); j++) {
      if ((spsf.componentsHO[j].n == spsi.componentsHO[i].n) &&
          (spsf.componentsHO[j].l == spsi.componentsHO[i].l)) {
        double I = ChargeDistributions::GetRadialMEHO(
            spsf.componentsHO[j].n, spsf.componentsHO[j].l, 0,
            spsi.componentsHO[i].n, spsi.componentsHO[i].l, nu);
        double r2 = ChargeDistributions::GetRadialMEHO(
            spsf.componentsHO[j].n, spsf.componentsHO[j].l, 2,
            spsi.componentsHO[i].n, spsi.componentsHO[i].l, nu);

        result += std::pow(spsf.componentsHO[j].C * spsi.componentsHO[i].C, 2.) * I *
                  (I - 2. * epsilon * r2);
        C += std::pow(spsf.componentsHO[j].C * spsi.componentsHO[i].C, 2.);
      }
    }
  }
  result *= (1. + 6. / 5. * epsilon * R * R) / C;

  return result;*/
  return 1;
}

double BSG::SpectralFunctions::RelativisticCorrection(double W, double W0, int Z,
                                                 double R, int betaType,
                                                 NHL::BetaDecayType decayType) {
  if (decayType == NHL::BetaDecayType::FERMI) {
    double Wb = W + betaType * 3 * fine_structure_const * Z / (2. * R);
    double pb = std::sqrt(Wb * Wb - 1.);
    double H2 = -std::pow(pb * R, 2.) / 6.;
    double D1 = Wb * R / 3.;
    double D3 =
        -Wb * R * std::pow(pb * R, 2.) / 30 - betaType * fine_structure_const * Z / 10.;
    double d1 = R / 3.;
    double d3 = -R * std::pow(pb * R, 2.) / 30.;
    double N1 = (W0 - W) * R / 3.;
    double N2 = -std::pow((W0 - W) * R, 2.) / 6.;
    double N3 = std::pow((W0 - W) * R, 3.) / 30;

    double gamma = std::sqrt(1. - std::pow(fine_structure_const * Z, 2.));

    double Vf2, Vf3;
    double Af2, Af3;

    Vf2 = -2. * (D1 + N1) + 2 * gamma / W * d1;
    Vf3 = -2. * (D3 + N1 * H2 - N2 * D1 - N3) + 2 * gamma / W * (d3 - N2 * d1);

    Af2 = 2. / std::sqrt(3) * (D1 + N1) - 2. / std::sqrt(3) * gamma / W * d1;
    Af3 = 2. * std::sqrt(2. / 3.) * (D1 - N1) -
          2. * std::sqrt(2. / 3.) * gamma / W * d1;

    double mismatch = W0 - betaType * 2.5 + betaType * 6. / 5. * fine_structure_const * Z / R;

    //TODO Check why Af are unused
    return 1. - 3. / 10 * R * mismatch * Vf2 - 3. / 28. * R * mismatch * Vf3;
  } else {
    return 1.;
  }
}

struct my_L0_params {double a; double b; double W; int Z; int betaType; std::array<double, 7> aPos; std::array<double, 7> aNeg;};

double L0Integral(double r, void * p) {
  struct my_L0_params * params = (struct my_L0_params *)p;
  double a = (params->a);
  double b = (params->b);
  double W = (params->W);
  int Z = (params->Z);
  int betaType = (params->betaType);
  std::array<double, 7> aPos = (params->aPos);
  std::array<double, 7> aNeg = (params->aNeg);
  double gamma = std::sqrt(1. - std::pow(fine_structure_const * Z, 2.));

  return std::pow(r, 3.) *
   abs(NHL::GetDerivDeformedChargeDist(r, a, b))
   * BSG::SpectralFunctions::L0Correction(W, Z, r, betaType, aPos, aNeg)
   * std::pow(r, 2. * (gamma - 1));
}

double BSG::SpectralFunctions::DeformationCorrection(double W, double W0, int Z,
                                                double R, double beta2,
                                                int betaType, std::array<double, 7> aPos,
                                                std::array<double, 7> aNeg) {
  double bOverA = NHL::CalcBoverA(beta2);
  double a, b;

  a = R * std::sqrt(3. / (bOverA * bOverA + 2.));
  b = bOverA * a;

  double V0s = 3. * fine_structure_const * Z / 2. / R;
  double V0d = 0.;

  if (beta2 > 0) {
    V0d = 3. * fine_structure_const * Z / 2. / a / std::sqrt(bOverA * bOverA - 1.) *
          std::log(std::sqrt(bOverA * bOverA - 1) + bOverA);
  } else if (beta2 < 0) {
    V0d = 3. * fine_structure_const * Z / 2. / a / std::sqrt(1. - bOverA * bOverA) *
          asin(std::sqrt(1 - bOverA * bOverA));
  }
  double etaS =
      1. / 6. * (std::pow(W0 - W, 2.) + std::pow(W + betaType * V0s, 2.) - 1.);
  double etaD =
      1. / 6. * (std::pow(W0 - W, 2.) + std::pow(W + betaType * V0d, 2.) - 1.);
  double DC0 = (1 - etaD * 3. / 5. * R * R) / (1 - etaS * 3. / 5. * R * R);

  double gamma = std::sqrt(1. - std::pow(fine_structure_const * Z, 2.));
  double prefact = 4. / 3. * pi * std::pow(R, 2. * (1 - gamma));

  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  double result, error;

  struct my_L0_params params = {a, b, W, Z, betaType, aPos, aNeg};

  gsl_function F;
  F.function = &L0Integral;
  F.params = &params;

  gsl_integration_qags (&F, a, b, a, 1e-7, 1000,
                        w, &result, &error);

  double DFS = prefact * result / L0Correction(W, Z, R, betaType, aPos, aNeg);

  gsl_integration_workspace_free (w);

  return DC0 * DFS;
}

double BSG::SpectralFunctions::NeutrinoRadiativeCorrection(double Wv) {
  double h = 0.;
  double pv = std::sqrt(Wv * Wv - 1);
  double beta = pv / Wv;
  h += 3 * std::log(proton_mass_c2 / electron_mass_c2) + 23 / 4. +
       8 / beta * Spence(2 * beta / (1 + beta)) +
       8 * (std::atan(beta) / beta - 1) * std::log(2 * Wv * beta) +
       4 * std::atan(beta) / beta *
           ((7 + 3 * beta * beta) / 8 - 2 * std::atan(beta));
  return 1 + fine_structure_const / 2 / pi * h;
}

double BSG::SpectralFunctions::Spence(double x) { return -gsl_sf_dilog(x); }

double BSG::SpectralFunctions::AtomicScreeningCorrection(double W, int Z,
                                                    int betaType) {
  std::vector<double> Aby, Bby;

  screening::PotParam(Z - 1 * betaType, Aby, Bby);

  double l = 2 * (Aby[0] * Bby[0] + Aby[1] * Bby[1] + Aby[2] * Bby[2]);

  double p = std::sqrt(W * W - 1);

  double Wt = W - betaType * 0.5 * fine_structure_const * (Z - betaType) * l;

  std::complex<double> pt;

  pt = 0.5 * p +
       0.5 * std::sqrt(std::complex<double>(p * p -
                                            betaType * 2 * fine_structure_const * Z * Wt * l));

  double y = betaType * fine_structure_const * Z * W / p;
  std::complex<double> yt = betaType * fine_structure_const * Z * Wt / pt;

  double gamma = std::sqrt(1. - pow(fine_structure_const * Z, 2.));

  /*cout << "V: " << W-Wt << endl;
  cout << W << " " << yt.real() << " " << yt.imag() << endl;
  cout << W << " " << pt.real() << " " << pt.imag() << endl;*/

  gsl_sf_result magn;
  gsl_sf_result phase;
  gsl_sf_lngamma_complex_e(gamma, y, &magn, &phase);

  gsl_sf_result magnT;
  gsl_sf_result phaseT;
  gsl_sf_lngamma_complex_e(gamma - yt.imag(), yt.real(), &magnT, &phaseT);
  // now we have what we wAt in magn.val

  double first = Wt / W;
  double second = std::exp(2 * (magnT.val - magn.val));

  gsl_sf_lngamma_complex_e(gamma - 2 * pt.imag() / l, 2 * pt.real() / l, &magnT,
                           &phaseT);
  gsl_sf_lngamma_complex_e(1, 2 * p / l, &magn, &phase);
  double third = std::exp(2 * (magnT.val - magn.val));
  double fourth = std::exp(-pi * y);
  double fifth = std::pow(2 * p / l, 2 * (1 - gamma));

  return first * second * third * fourth * fifth;
}

double BSG::SpectralFunctions::AtomicExchangeCorrection(double W, std::array<double, 9> exPars) {
  double E = W - 1;

  return 1 + exPars[0] / E + exPars[1] / E / E +
         exPars[2] * std::exp(-exPars[3] * E) +
         exPars[4] * sin(std::pow(W - exPars[6], exPars[5]) + exPars[7]) /
             std::pow(W, exPars[8]);
}

double BSG::SpectralFunctions::AtomicMismatchCorrection(double W, double W0, int Z,
                                                   int A, int betaType) {
  double dBdZ2 = (44.200 * std::pow(Z - betaType, 0.41) +
                  2.3196E-7 * std::pow(Z - betaType, 4.45)) * eV /
                 electron_mass_c2;

  double K = -0.872 + 1.270 * std::pow(Z, 0.097) + 9.062E-11 * std::pow(Z, 4.5);
  double vp = std::sqrt(1 - 1 / W / W);
  double l = 1.83E-3 * K * Z / vp;
  double M = A * NHL::nucleon_mass_bu;
  // assume as A average the recoil velocity at half-momentum trAsfer ~
  // std::sqrt(W0^2-1)/2
  double vR = std::sqrt(1 - M * M / (M * M + (W0 * W0 - 1) / 4.));

  double psi2 = 1 + 2 * fine_structure_const / vp * (std::atan(1 / l) - l / 2 / (1 + l * l));

  double C0 = -fine_structure_const * fine_structure_const * Z * fine_structure_const / vp * l / (1 + l * l) / psi2;

  double C1 = 2 * fine_structure_const * fine_structure_const * Z * vR / vp *
              ((0.5 + l * l) / (1 + l * l) - l * std::atan(1 / l)) / psi2;

  return 1 - 2 / (W0 - W) * (0.5 * dBdZ2 + 2 * (C0 + C1));
}
