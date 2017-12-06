#ifndef NILSSON_ORBITS
#define NILSSON_ORBITS

#include "Utilities.h"
#include "Constants.h"
#include "NuclearUtilities.h"
#include "spdlog.h"

#include <iostream>
#include <algorithm>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#define NDIM1 56
#define NDIM3 35
#define NDIM4 40

namespace NuclearStructure {
/**
 * This ia a C++ port of the original Fortran IV code from 
 * Computer Physics Communications 6 (1973) 30
 * by Hird et al.
 * 
 * Here follows the original documentation
 * 
 * ABOVNILSSON ORBITS. NILSSON ORBITS FOR A PARTICLE IN A WOODS-SAXON
 * 1   POTENTIAL WITH Y20 AND Y40 DEFORMATIONS, AND COUPLED TO CORE
 * 2   ROTATIONAL STATES.  HIRD, B.
 * REF. IN COMP. PHYS. COMMUN. 6 (1973) 30
 * NILSSON ORBITS
 * FOR A PARTICLE IN A WOODS-SAXON POTENTIAL WITH Y(2,0) AND Y(4,0)
 * DEFORMATIONS, AND COUPLED TO CORE ROTATIONAL STATES.

 * METHOD - MATRIX DIAGONALIZATION.
 *    FIRST DIAGONALIZATION IS OF SPHERICAL WOODS-SAXON HAMILTONIAN
 *         WITH HARMONIC OSCILLATOR BASIS.
 *    SECOND DIAGONALIZATION IS OF DEFORMED HAMILTONIAN WITH
 *         SPHERICAL WOODS-SAXON EIGENSTATES AS BASIS.
 *    THIRD DIAGONALIZATION IS OF BAND MIXING AND CORE HAMILTONIAN
 *         WITH DEFORMED PARTICLE EIGENSTATES AS A BASIS.
 *
 * REFERENCES:-
 * S.G.NILSSON, 1955 MAT.FYS.MEDD.DAN.VID.SELSK. 29 NO 16.
 * A.K.KERMAN, 1956 MAT.FYS.MEDD.DAN.VID.SELSK. 30 NO 15.
 * K.T.HECHT AND G.R.SATCHLER, NUCLEAR PHYSICS, 32(1962)286
 * W.SCHOLZ AND F.B.MALIK, PHYS.REV., 147(1966)836.

 * NOTATION:-
 * ENERGIES IN MEV.  LENGTHS IN FM. MASS IN A.M.U. Z IN ELECTRONIC
 * CHARGE UNITS.
 * V0,R0,A0,V0S = WOODS-SAXON PARAMETERS.
 * AMU = ATOMIC WEIGHT OF ODD PARTICLE + CORE.
 * Z = ZERO FOR A NEUTRON STATE; = CORE CHARGE FOR A PROTON STATE.
 * NMAX = MAXIMUM RADIAL QUANTUM NUMBER IN HARMONIC OSCILLATOR BASIS.
 * NMAX = ODD;EVEN FOR NEGATIVE;POSITIVE PARITY     NMAX < 14
 * BETA2 = SPHEROIDAL DEFORMATION
 * BETA4 = HEXADECAPOLE DEFORMATION.
 * ROT = UNIT OF ROTATIONAL ENERGY OF THE CORE.
 * ROTI(I) = ROTATIONAL ENERGY FOR THIS VALUE OF I ONLY.
 * ROTO(OMEGA) = ROTATIONAL ENERGY FOR THIS OMEGA ONLY.
 * SPIN = MAXIMUM VALUE OF TOTAL SPIN I OF SYSTEM.
 * ALL EIGENSTATES WITH SPINS UP TO INPUT VALUE OF SPIN ARE FOUND.
 * N = RADIAL QUANTUM NUMBER.
 * L = ANGULAR MOMENTUM QUANTUM NUMBER.
 * LA = LAMBDA = COMPONENT OF L ALONG SYMMETRY AXIS.
 * IX2 = 2*SIGMA = 2*COMPONENT OF INTRINSIC SPIN ALONG SYMMETRY AXIS.
 * JX2 = 2*J, WHERE J = L + INTRINSIC SPIN.
 * SW(1,NN) = TABLE OF <N',L!F(R)!N,L>
 * SW(2,NN) = TABLE OF <N',L!(1/R)*D(F(R))/D(R)!N,L>.
 * SDW(NNP) = TABLE OF <N',L'!D(F(R)/D(R)!N,L>.                      
 */
namespace nilsson {

using std::cout;
using std::endl;

/**
 * Value of harmonic oscillator function at radius x
 * 
 * @param n principal quantum number
 * @param l orbital quantum number
 * @param x radius in atomic units
 */
inline double V(int n, int l, double x) {
  int m = (n - l) / 2;
  double V = 1.0;
  double F = 1.0;
  int mK = 1;
  int mM = 1;
  if (m > 0) {
    for (int k = 1; k <= m; k++) {
      mK *= m + 1 - k;
      F *= -2. * x / k;
      mM *= 2 * n - 4 * m + 1 + 2 * k;
      V += F * mK / mM;
    }
  }
  return V;
}

/**
 * Calculate the normalization constant for a harmonic oscillator wave function
 * 
 * @param n principal quantum number
 * @param l orbital quantum number
 */
inline double VNORM(int n, int l) {
  int minus = (n - l) / 2;
  double eMinus = 1 + 4 * (minus / 2) - 2 * minus;
  double vNorm = 2. / std::pow(M_PI, 0.25) * eMinus;
  int k = (n + 1) / 2 - n / 2 + 1;
  if (k > 1) {
    vNorm *= std::sqrt(2. / 3.);
  }
  if (n > 1) {
    for (int j = k; j <= n; j += 2) {
      vNorm *= std::sqrt((j + k + 1.) / (j - k + 2.));
    }
  }
  if (l > 1) {
    for (int j = k; j <= l; j += 2) {
      vNorm *= 2. * std::sqrt((j + n + 2.) * (n - j + 1.)) /
               ((2. * j + 3.) * (2. * j + 1.));
    }
  }
  return vNorm;
}

/**
 * Calculate the radial integrals for all harmonic oscillator functions in a Woods-Saxon potential
 * 
 * @param V0 depth of the Woods-Saxon potential
 * @param R nuclear radius
 * @param A0 
 * @param V0S strength of the pion exchange term
 * @param A mass number
 * @param Z proton number
 * @param nMax maximum number of oscillator shells
 * @param SW array containing values for radial integrals in Woods-Saxon potential
 * @param SDW array containing values for radial integrals in deformed Woods-Saxon potential
 */
inline void WoodsSaxon(double V0, double R, double A0, double V0S, double A,
                       double Z, int nMax, double SW[2][84], double SDW[462]) {
  auto dbl = spdlog::get("debug_file");
  dbl->debug("Entered WoodsSaxon");
  double FINT[3] = {};
  double S[3][4] = {};
  int II = 0;
  int JJ = 0;
  int NDX = 300;
  int N = NDX / 4 - 2;
  double AO = A0;
  // Compton wavelength of the pion=  1.4133 fm
  double pionComp = HBAR / PION_MASS_KEV / SPEED_OF_LIGHT;
  pionComp = 1.4133;
  double VOS = V0S * pionComp * pionComp;
  double TWONU = 0.154666 * std::sqrt(V0) / R;
  double T = TWONU * 20.899 * A / (A - 1.);
  // e^2 = 1.439978 MeV*fm
  double ZZ = Z * 1.439978;
  double ZCONST = 1.5 * ZZ / R;
  double ZIN = 0.5 * ZZ / std::pow(R, 3.);
  double DX = (R + 8. * AO) / NDX;
  double DFO = std::exp(DX / AO);
  int nMin = nMax % 2;

  dbl->debug("Radial integrals for Woods-Saxon potential");
  for (int NI = nMin; NI <= nMax; NI += 2) {
    for (int NJ = nMin; NJ <= NI; NJ += 2) {
      for (int LI = nMin; LI <= NI; LI += 2) {
        for (int LJ = nMin; LJ <= NJ; LJ += 2) {
          int KK = 0;
          if (LI != LJ) {
            KK = 2;
          }
          double r = DX;
          double FO = DFO / std::exp(R / AO);
          for (int i = 0; i < 3; i++) {
            FINT[i] = 0.0;
          }
          int M = -N - 2;
          while (true) {
            double r2 = r * r * TWONU;
            double PSI2 = V(NI, LI, r2) * V(NJ, LJ, r2) *
                          std::pow(r2, (LI + LJ) / 2) * r * r / std::exp(r2);
            double F = 1. / (1. + FO);
            double DFDR = -FO * F * F / AO;
            for (int i = 0; i < 4; i++) {
              if (r <= R) {
                S[0][i] = -PSI2 * (V0 * F + T * r2 + ZIN * r * r - ZCONST);
              } else {
                S[0][i] = -PSI2 * (V0 * F + T * r2 - ZZ / r);
              }
              S[1][i] = PSI2 * DFDR / r;
              S[2][i] = PSI2 * DFDR;
              FO *= DFO;
              r += DX;
              for (int j = KK; j < 3; j++) {
                FINT[j] += S[j][i];
              }
            }
            M += 2;
            if (M / N < 0) {
              for (int i = 1; i < 4; i++) {
                for (int l = i; l < 4; l++) {
                  int k = 3 + i - l;
                  for (int j = KK; j < 3; j++) {
                    S[j][k] -= S[j][k - 1];
                  }
                }
              }
              for (int j = KK; j < 3; j++) {
                FINT[j] += -S[j][0] / 2. + S[j][1] / 12. - S[j][2] / 24. +
                           19. * S[j][3] / 720.;
              }
              if (DX <= 0) {
                break;
              }
            } else if (M / N > 0) {
              r += 4. * DX;
              DX = -DX;
              FO *= std::pow(DFO, 4.);
              DFO = 1. / DFO;
              M = -N - 2;
            }
          }
          DX = -DX;
          DFO = 1. / DFO;
          double X = VNORM(NI, LI) * VNORM(NJ, LJ);
          if (KK != 2) {
            SW[0][II] = FINT[0] * X * DX * std::pow(TWONU, 1.5);
            if (NI == NJ) {
              // Harmonic oscillator energy
              SW[0][II] += 2. * T * (NI + 1.5);
            }
            SW[1][II] = FINT[1] * X * DX * std::pow(TWONU, 1.5) * VOS;
            II++;
          }
          SDW[JJ] = FINT[2] * X * DX * std::pow(TWONU, 1.5) * V0 * R;
          JJ++;
          if (KK == 0) {
            dbl->debug("< {},{} | {},{} >", NI, LI, NJ, LJ);
            dbl->debug("    {} {} =    {}\t {}\t {}", II, JJ, SW[0][II - 1], SW[1][II - 1], SDW[JJ - 1]);
          } else if (KK == 2) {
            dbl->debug("< {},{} | {},{} >", NI, LI, NJ, LJ);
            dbl->debug("      {} = \t\t\t\t\t{}", JJ, SDW[JJ - 1]);
          }
        }
      }
    }
  }
  dbl->debug("Leaving WoodsSaxon");
}

// Calculate eigenvalues & vectors for a real symmetric FORTRAN matrix A, only
// upper half of A is used
/**
 * Calculate the eigen values and eigen vectors of a matrix
 * only the upper half of A is used by default
 * 
 * @param A pointer to an array containing the matrix elements in symmetric FORTRAN style
 * @param dim dimension of the matrix
 * @param eVecs pointer to an array in which to place the eigenvectors
 * @param eVals reference to a vector in which to put the eigenvalues
 * @param onlyUpper boolean to say whether only the upper part was given
 */
inline void Eigen(double* A, int dim, double* eVecs, std::vector<double>& eVals,
                  bool onlyUpper = true) {
  auto dbl = spdlog::get("debug_file");
  dbl->debug("Entered Eigen");
  gsl_matrix* aNew = gsl_matrix_alloc(dim, dim);
  // Loop over upper half of matrix
  for (int j = 0; j < dim; j++) {
    for (int i = 0; i <= j; i++) {
      // FORTRAN matrices are stored column-major
      int index = dim * j + i;
      if (onlyUpper) {
        index = j * (j + 1) / 2 + i;
      }
      gsl_matrix_set(aNew, i, j, A[index]);
      if (i != j) {
        gsl_matrix_set(aNew, j, i, A[index]);
      }
    }
  }
  gsl_matrix* eVec = gsl_matrix_calloc(dim, dim);
  gsl_vector* eVal = gsl_vector_calloc(dim);
  int size = 4 * dim;
  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv(aNew, eVal, eVec, w);
  gsl_eigen_symmv_free(w);
  // Sort eigenvalues and eigenvalues according to descending eigenvalue
  for (int i = 0; i < dim; i++) {
    for (int j = i; j < dim; j++) {
      if (gsl_vector_get(eVal, j) > gsl_vector_get(eVal, i)) {
        double x = gsl_vector_get(eVal, i);
        gsl_vector_set(eVal, i, gsl_vector_get(eVal, j));
        gsl_vector_set(eVal, j, x);
        for (int k = 0; k < dim; k++) {
          x = gsl_matrix_get(eVec, k, i);
          gsl_matrix_set(eVec, k, i, gsl_matrix_get(eVec, k, j));
          gsl_matrix_set(eVec, k, j, x);
        }
        i = 0;
        j = 0;
      }
    }
  }
  for (int i = 0; i < dim; i++) {
    eVals.push_back(gsl_vector_get(eVal, i));
  }
  gsl_vector_free(eVal);
  for (int j = 0; j < dim; j++) {
    for (int i = 0; i < dim; i++) {
      // FORTRAN matrices are stored column-major
      eVecs[dim * j + i] = gsl_matrix_get(eVec, i, j);
    }
  }
  gsl_matrix_free(eVec);
  dbl->debug("Leaving Eigen");
}

/**
 * Calculate the single particle eigenstates of the Woods-Saxon potential
 * Will calculate the deformed case when deformation parameters are non-zero
 * 
 * @param spin nuclear spin
 * @param beta2 quadrupole deformation
 * @param beta4 hexadecupole deformation
 * @param beta6 beta6 deformation
 * @param V0 depth of the Woods-Saxon potential
 * @param R nuclear radius in atomic units
 * @param A0 
 * @param V0S strength of the pion-spin-exchange
 * @param A mass number
 * @param Z proton number
 * @param nMax maximum number of oscillator shells
 * @returns vector a SingleParticleState objects of all bound eigenstates in the potential
 */
inline std::vector<SingleParticleState> Calculate(
    double spin, double beta2, double beta4, double beta6, double V0, double R,
    double A0, double V0S, double A, double Z, int nMax) {
  auto dbl = spdlog::get("debug_file");
  dbl->debug("Entered Calculate");
  double SW[2][84] = {};
  double SDW[462] = {};
  int N[NDIM1] = {};
  int L[NDIM1] = {};
  // Matrix elements of Hamiltonian in chosen basis, real & symmetric
  double hamM[NDIM4 * (NDIM4 + 1) / 2] = {};
  double eVecs[NDIM4 * NDIM4] = {};
  double defExpCoef[NDIM4][NDIM1] = {};
  double B[NDIM4 * (NDIM4 + 1) / 2] = {};
  double D[NDIM4 * (NDIM4 + 1) / 2] = {};
  double F[NDIM4 * (NDIM4 + 1) / 2] = {};
  int LA[NDIM1] = {};
  int IX2[NDIM1] = {};
  int JX2[NDIM1] = {};
  double sphExpCoef[NDIM3][NDIM1] = {};
  double eValsWS[NDIM3] = {};
  int NDOM[NDIM3] = {};
  int NDOMK[NDIM3] = {};
  int LK[NDIM3] = {};
  int LKK[NDIM3] = {};
  int KN[7][NDIM3] = {};
  int K1[7] = {};
  double eValsDWS[NDIM4] = {};
  int K2[NDIM4] = {};
  int K3[NDIM4] = {};
  int K4[NDIM4] = {};

  std::vector<SingleParticleState> states;

  WoodsSaxon(V0, R, A0, V0S, A, Z, nMax, SW, SDW);

  dbl->debug("Past WoodsSaxon");

  int II = 0;
  int K = 0;

  int nMin = nMax % 2 + 1;
  int nMaxP1 = nMax + 1;
  // Diagonalize each L-value separately
  for (int LI = nMin; LI <= nMaxP1; LI += 2) {
    int NIM = II;
    int NI = II + 1;
    int NK = K + 1;

    // Set up quantum numbers of the harmonic oscillator basis
    for (int NN = LI; NN <= nMaxP1; NN += 2) {
      for (int I = 1; I <= 2; I++) {
        N[II] = NN - 1;
        L[II] = LI - 1;
        LA[II] = 2 - I;
        IX2[II] = 2 * I - 3;
        JX2[II] = 2 * L[II] + IX2[II];
        if (L[II] >= LA[II]) {
          II++;
        }
      }
    }
    if (II == NIM) {
      break;
    }

    // Set up matrix
    int KK = 0;
    for (int i = NI; i <= II; i++) {
      for (int j = NI; j <= i; j++) {
        hamM[KK] = 0.0;
        if (JX2[j - 1] == JX2[i - 1]) {
          int NN =
              (N[i - 1] / 2) * (N[i - 1] / 2 + 1) * (N[i - 1] / 2 + 2) / 6 +
              (N[j - 1] / 2) * (N[j - 1] / 2 + 1) / 2 + L[i - 1] / 2;
          hamM[KK] =
              SW[0][NN] + SW[1][NN] * IX2[j - 1] * (L[j - 1] + LA[j - 1]);
        }
        KK++;
      }
    }
    std::vector<double> eVals;
    Eigen(hamM, II - NIM, eVecs, eVals);

    int index = 0.0;
    for (int i = NI; i <= II; i++) {
      double eigenValue = eVals[index];
      if (eigenValue <= 10.0) {
        if (K >= NDIM3 - 1) {
          dbl->warn("Dimensioned space inadequate");
          return states;
        }
        eValsWS[K] = eigenValue;
        LK[K] = L[i - 1];
        K3[K] = JX2[i - 1];
        for (int J = NI; J <= II; J++) {
          int N0 = (II - NIM) * (i - NI) + J - NIM;
          sphExpCoef[K][J - 1] = eVecs[N0 - 1];
        }
        NDOM[K] = 0;
        double max = 0.0;
        for (int J = NI; J <= II; J++) {
          if (std::abs(sphExpCoef[K][J - 1]) > max) {
            NDOM[K] = N[J - 1];
            max = std::abs(sphExpCoef[K][J - 1]);
          }
        }
        K++;
      }
      index++;
    }
    // II = size of harmonic oscillator basis
    // K = size of Woods-Saxon basis for deformed state diagonalization
    if (NK == K + 1) {
      II = NIM;
    }
  }
  if (K == 0) {
    dbl->warn("No harmonic oscillator single particle states below 10 MeV.");
    return states;
  }

    //spdlog::get("nme_result_file")->info("Spherical Woods-Saxon expansion in Harmonic Oscillator basis.");
    for (int KKK = 1; KKK <= K; KKK += 12) {
      int KKKK = std::min(K, KKK + 11);
      cout << "Energy (MeV): ";
      for (int i = KKK; i <= KKKK; i++) {
        cout << eValsWS[i - 1] << " ";
      }
      cout << endl;
      cout << "Spin: \t\t";
      for (int i = KKK; i <= KKKK; i++) {
        cout << K3[i - 1] << "/2 \t";
      }
      cout << endl;
      cout << "Coefficients: " << endl;
      for (int J = 1; J <= II; J++) {
        cout << "| " << N[J - 1] << " " << JX2[J - 1] << "/2 > ";
        for (int L = KKK; L <= KKKK; L++) {
          cout << sphExpCoef[L - 1][J - 1] << "\t\t ";
        }
        cout << endl;
      }
      cout << "Dominant N: ";
      for (int i = KKK; i <= KKKK; i++) {
        cout << NDOM[i - 1] << "\t";
      }
      cout << endl;
    }
  dbl->debug("Past spherical case");
  if (!(beta2 == 0 && beta4 == 0 && beta6 == 0)) {
    int ISX2 = std::abs(2 * spin);
    int ISPIN = (ISX2 + 1) / 2;

    int KK = 0;
    int KKKK = 0;
    int IOM = 1;
    for (int IIOM = 1; IIOM <= ISPIN; IIOM++) {
      IOM = -IOM - 2 * (int)(std::pow(-1., IIOM));
      int IIIOM = 4 * (std::abs(IOM) / 4) + 2 - nMin;
      int KKK = 0;
      // Loop over all spherical states with E < 10.0 MeV
      for (int I = 0; I < K; I++) {
        LKK[KKK] = LK[I];
        KN[IIOM - 1][KKK] = I + 1;
        NDOMK[KKK] = NDOM[I];
        double X = 0.0;
        // Loop over the full spherical basis
        for (int J = 0; J < II; J++) {
          if (L[J] == LKK[KKK]) {
            double cg = utilities::ClebschGordan(
                2 * L[J], 1, JX2[J], 2 * (LA[J] + IIOM - 1), IX2[J],
                2 * (LA[J] + IIOM - 1) + IX2[J]);
            double X1 = sphExpCoef[I][J] * cg;
            cg = utilities::ClebschGordan(2 * L[J], 1, JX2[J] - 2 * IX2[J],
                                          2 * (LA[J] + IIOM - 1), IX2[J],
                                          2 * (LA[J] + IIOM - 1) + IX2[J]);
            if (JX2[J] > IIIOM) {
              X1 += sphExpCoef[I][J - IX2[J]] * cg;
            }
            defExpCoef[KKK][J] = X1;
            X += X1 * X1;
          }
        }
        if (X < 0.1) {
          KKK--;
        }
        KKK++;
      }
      KKKK += KKK;
      if (KKKK > NDIM4) {
        cout << "Problem KKKK > NDIM4" << endl;
        return states;
      }
      K1[IIOM - 1] = KKK;
      if (KKK == 0) {
        break;
      }
      for (int I = 0; I < KKK; I++) {
        for (int J = 0; J <= I; J++) {
          B[KK] = 0.0;
          D[KK] = 0.0;
          F[KK] = 0.0;
          for (int N1 = 0; N1 < II; N1++) {
            for (int N2 = 0; N2 < II; N2++) {
              if (!(L[N1] != LKK[I] || L[N2] != LKK[J] ||
                    (L[N1] - L[N2]) / 5 != 0 || LA[N1] != LA[N2])) {
                int NI = N[N1] / 2 + 1;
                int NJ = N[N2] / 2 + 1;
                int LI = L[N1] / 2 + 1;
                int LJ = L[N2] / 2 + 1;
                if (NI < NJ) {
                  int NN = NI;
                  NI = NJ;
                  NJ = NN;
                  NN = LI;
                  LI = LJ;
                  LJ = NN;
                }
                int NNP = NI * (NI - 1);
                NNP = (3 * NNP * NNP + 2 * NNP * (2 * NI - 1)) / 24 +
                      NI * NJ * (NJ - 1) / 2 + (LI - 1) * NJ + LJ;
                double X = defExpCoef[I][N1] * SDW[NNP - 1] * defExpCoef[J][N2];
                int LP = L[N1];
                int LAP = LA[N1] + IIOM - 1;
                int LL = L[N2];
                int LLA = LA[N2] + IIOM - 1;
                B[KK] +=
                    X * utilities::SphericalHarmonicME(LP, LAP, 2, 0, LL, LLA);
                D[KK] +=
                    X * utilities::SphericalHarmonicME(LP, LAP, 4, 0, LL, LLA);
                F[KK] +=
                    X * utilities::SphericalHarmonicME(LP, LAP, 6, 0, LL, LLA);
              }
            }
          }
          KK++;
        }
      }
    }

    K = 0;
    KK = 0;
    IOM = 1;
    for (int IIOM = 1; IIOM <= ISPIN; IIOM++) {
      IOM = -IOM - 2 * (int)(std::pow(-1., IIOM));
      int KKK = K1[IIOM - 1];
      if (KKK != 0) {
        int NK = 0;
        for (int I = 1; I <= KKK; I++) {
          for (int J = 1; J <= I; J++) {
            NK++;
            hamM[NK - 1] = beta2 * B[KK] + beta4 * D[KK] + beta6 * F[KK];
            KK++;
          }
          int N0 = KN[IIOM - 1][I - 1];
          hamM[NK - 1] += eValsWS[N0 - 1];
        }
        std::vector<double> eVals;
        Eigen(hamM, KKK, eVecs, eVals);
        int N0 = 0;
        for (int MU = 1; MU <= KKK; MU++) {
          N0 += MU;
          K2[K] = IOM;
          K3[K] = std::abs(IOM);
          K4[K] = KKK - MU + 1;
          eValsDWS[K] = eVals[MU - 1];
          for (int J = 0; J < II; J++) {
            defExpCoef[K][J] = 0.0;
            NK = KKK * (MU - 1);
            for (int NU = 0; NU < KKK; NU++) {
              int I = KN[IIOM - 1][NU];
              NK++;
              defExpCoef[K][J] += eVecs[NK - 1] * sphExpCoef[I - 1][J];
            }
          }
          NDOMK[K] = 0;
          double max = 0.0;
          for (int j = 0; j < II; j++) {
            if (std::abs(defExpCoef[K][j]) > max) {
              NDOMK[K] = N[j];
              max = std::abs(defExpCoef[K][j]);
            }
          }
          K++;
        }
      }
    }
      cout << "Deformed states: No band mixing   Beta2: " << beta2
           << " Beta4: " << beta4 << " Beta6: " << beta6 << endl;
      for (int KKK = 1; KKK <= K; KKK += 12) {
        int KKKK = std::min(K, KKK + 11);
        cout << "Energy (MeV): ";
        for (int i = KKK; i <= KKKK; i++) {
          cout << eValsDWS[i - 1] << " ";
        }
        cout << endl;
        cout << "*OMEGA |MU> ";
        for (int i = KKK; i <= KKKK; i++) {
          cout << K3[i - 1] << "/2|" << K4[i - 1] << "> \t";
        }
        cout << endl;
        for (int j = 1; j <= II; j++) {
          cout << "| " << N[j - 1] << ", " << JX2[j - 1] << "/2 > " << LA[j - 1]
               << " ";
          for (int l = KKK; l <= KKKK; l++) {
            cout << defExpCoef[l - 1][j - 1] << "\t\t";
          }
          cout << endl;
        }
        cout << "Dominant N: ";
        for (int i = KKK; i <= KKKK; i++) {
          cout << NDOMK[i - 1] << "\t";
        }
        cout << endl;
      }
  }  // end of deformation only part
  for (int i = 0; i < K; i++) {
    SingleParticleState sps;
    if (beta2 == 0 && beta4 == 0) {
      sps.energy = eValsWS[i];
    } else {
      sps.energy = eValsDWS[i];
    }
    sps.dO = K3[i];
    sps.dK = K3[i];
    sps.nDom = NDOMK[i];
    int nZ = NDOMK[i] - (K3[i] - 1) / 2;
    for (int j = 0; j < K; j++) {
      if (eValsDWS[j] < sps.energy && NDOMK[j] == NDOMK[i] && K3[j] == K3[i]) {
        nZ--;
      }
    }
    sps.nZ = nZ;
    if ((NDOMK[i] - nZ) % 2 == 0) {
      if ((K3[i] + 1) / 2 % 2 == 0) {
        sps.lambda = (K3[i] + 1) / 2;
      } else {
        sps.lambda = (K3[i] - 1) / 2;
      }
    } else {
      if ((K3[i] + 1) / 2 % 2 == 0) {
        sps.lambda = (K3[i] - 1) / 2;
      } else {
        sps.lambda = (K3[i] + 1) / 2;
      }
    }
    sps.parity = 1 - 2 * (nMax % 2);
    for (int j = 0; j < II; j++) {
      WFComp wfc;
      // spectroscopic quantum number
      wfc.n = (N[j] - L[j]) / 2 + 1;
      wfc.l = L[j];
      wfc.s = (JX2[j] - L[j] * 2);
      if (beta2 == 0 && beta4 == 0) {
        wfc.C = sphExpCoef[i][j];
      } else {
        wfc.C = defExpCoef[i][j];
      }
      sps.componentsHO.push_back(wfc);
    }
    states.push_back(sps);
  }
  return states;
}

/**
 * Custom sort function to compare to SingleParticleState objects.
 * Sorts on energy of the state, checks whether lhs < rhs
 * 
 * @param lhs left-hand side
 * @param rhs right-hand side
 */
inline bool StateSorter(SingleParticleState const& lhs,
                        SingleParticleState const& rhs) {
  return lhs.energy < rhs.energy;
}

/**
 * Get a sorted vector of all bound Single particle states, even and odd parity
 * 
 * @param Z proton number
 * @param N neutron number
 * @param A mass number
 * @param dJ double of nuclear spin
 * @param R nuclear radius in atomic units
 * @param beta2 quadrupole deformation
 * @param beta4 hexadecupole deformation
 * @param beta6 beta6 deformation
 * @param V0 depth of the Woods-Saxon potential
 * @param A0
 * @param VS strength of the pion-exchange
 * @returns vector of all bound single particle states, sorted for increasing energy
 */
inline std::vector<SingleParticleState> GetAllSingleParticleStates(
    int Z, int N, int A, int dJ, double R, double beta2, double beta4,
    double beta6, double V0, double A0, double VS) {
  std::vector<SingleParticleState> evenStates =
      Calculate(6.5, beta2, beta4, beta6, V0, R, A0, VS, A, Z, 12);
  std::vector<SingleParticleState> oddStates =
      Calculate(6.5, beta2, beta4, beta6, V0, R, A0, VS, A, Z, 13);

  // Join all states
  std::vector<SingleParticleState> allStates;
  allStates.reserve(evenStates.size() + oddStates.size());
  allStates.insert(allStates.end(), evenStates.begin(), evenStates.end());
  allStates.insert(allStates.end(), oddStates.begin(), oddStates.end());

  // Sort all states according to energy
  std::sort(allStates.begin(), allStates.end(), &StateSorter);

  return allStates;
}

/**
 * Search for the deformation single particle state corresponding to a certain spin
 * within a certain energy threshold of the calculated ones
 * 
 * @param Z proton number
 * @param N neutron number
 * @param A mass number
 * @param dJ double of nuclear spin
 * @param R nuclear radius in atomic units
 * @param beta2 quadrupole deformation
 * @param beta4 hexadecupole deformation
 * @param beta6 beta6 deformation
 * @param V0 depth of the Woods-Saxon potential
 * @param A0
 * @param VS strength of the pion-exchange
 * @param dJreq double of the required spin
 * @param threshold maximum energy difference between the calculated state with
 *     the correct spin and that proposed as the one at the Fermi surface
 * @returns SingleParticleState object. If the correct state is not found, it returns
 *     the first SingleParticleState
 */
inline SingleParticleState CalculateDeformedSPState(int Z, int N, int A, int dJ,
                                                    double R, double beta2,
                                                    double beta4, double beta6,
                                                    double V0, double A0,
                                                    double VS, int dJreq,
                                                    double threshold) {
  std::vector<SingleParticleState> allStates = GetAllSingleParticleStates(
      Z, N, A, dJ, R, beta2, beta4, beta6, V0, A0, VS);

  int index = 0;
  if (beta2 == 0 && beta4 == 0 && beta6 == 0) {
    int nrParticles = Z + N;
    for (int i = 0; i < allStates.size(); i++) {
      if (nrParticles - (allStates[index].dO + 1) > 0) {
        nrParticles -= allStates[index].dO + 1;
        index++;
      }
    }
  } else {
    index = (Z + N - 1) / 2;
    if (threshold > 0) {
      double refEnergy = allStates[index].energy;
      index = 0;
      for (int i = 0; i < allStates.size(); i++) {
        if (std::abs(refEnergy - allStates[i].energy) <= threshold) {
          if (allStates[i].dO * allStates[i].parity == dJreq &&
              std::abs(allStates[i].energy - refEnergy) <
                  std::abs(allStates[index].energy - refEnergy)) {
            index = i;
          }
        }
      }
      if (index == 0) {
        cout
            << "ERROR: Couldn't find a correct spin state within the threshold."
            << endl;
      }
    }
  }

  cout << "Found single particle state. Energy: " << allStates[index].energy
       << endl;
  cout << "Spin: " << allStates[index].parity* allStates[index].dO << "/2 ["
       << allStates[index].nDom << allStates[index].nZ
       << allStates[index].lambda << "]" << endl;
  cout << "Orbital\tC" << endl;
  for (int j = 0; j < allStates[index].componentsHO.size(); j++) {
    cout << allStates[index].componentsHO[j].n
         << utilities::spectroNames[allStates[index].componentsHO[j].l]
         << allStates[index].componentsHO[j].l * 2 +
                allStates[index].componentsHO[j].s << "/2\t"
         << allStates[index].componentsHO[j].C << endl;
  }

  return allStates[index];
}

}
}
#endif
