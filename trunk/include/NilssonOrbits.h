#ifndef NILSSON_ORBITS
#define NILSSON_ORBITS

#include "Utilities.h"
#include "Constants.h"

#include <iostream>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#define NDIM1 56
#define NDIM3 35
#define NDIM4 40

namespace nilsson {

using std::cout;
using std::endl;

// Here dO stands for double Omega
// s is +-1 depending of whether it is j=l+-1/2
struct WFComp {
  double C;
  int l, s, dO;
};

struct NuclearState {
  double O, K;
  int parity;
  std::vector<WFComp> states;
};

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

inline void WoodsSaxon(double V0, double R, double A0, double V0S, double A,
                       double Z, int nMax, double SW[2][84], double SDW[462]) {
  cout << "In Woods-Saxon" << endl;
  double FINT[3] = {};
  double S[3][4] = {};
  int II = 0;
  int JJ = 0;
  int NDX = 300;
  int N = NDX / 4 - 2;
  double AO = A0;
  // Compton wavelength of the pion=  1.4133 fm
  double pionComp = hbar / pionMasskeV / speedOfLight;
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
            }
            else if (M / N > 0) {
              r += 4. * DX;
              DX = -DX;
              FO *= std::pow(DFO, 4.);
              DFO = 1. / DFO;
              M = -N - 2;
            }
          }
          DX = -DX;
          DFO = 1./DFO;
          double X = VNORM(NI, LI) * VNORM(NJ, LJ);
          if (KK != 2) {
            SW[0][II] = FINT[0] * X * DX * std::pow(TWONU, 1.5);
            if (NI == NJ) {
              //Harmonic oscillator energy
              SW[0][II] += 2. * T * (NI + 1.5);
            }
            SW[1][II] = FINT[1] * X * DX * std::pow(TWONU, 1.5) * VOS;
            II++;
          }
          SDW[JJ] = FINT[2] * X * DX * std::pow(TWONU, 1.5) * V0 * R;
          JJ++;
          if (KK == 0) {
            cout << "< " << NI << "," << LI << " | " << NJ << "," << LJ << " >";
            cout << "    " << II << " " << JJ << " =    " << SW[0][II - 1]
                 << " " << SW[1][II - 1] << " " << SDW[JJ - 1];
            cout << endl;
          } else if (KK == 2) {
            cout << "< " << NI << "," << LI << " | " << NJ << "," << LJ << " >";
            cout << "      " << JJ << " =    " << SDW[JJ - 1];
            cout << endl;
          }
        }
      }
    }
  }
}

// Calculate eigenvalues & vectors for a real symmetric FORTRAN matrix A
inline void Eigen(double* A, int rank, double* eVecs,
                  std::vector<double>* eVals) {
  gsl_matrix* aNew = gsl_matrix_alloc(NDIM4, NDIM4);
  for (int i = 0; i < NDIM4; i++) {
    for (int j = i; j < NDIM4; j++) {
      // FORTRAN matrices are stored column-major
      gsl_matrix_set(aNew, i, j, A[NDIM4 * i + j]);
      if (i != j) {
        gsl_matrix_set(aNew, j, i, A[NDIM4 * i + j]);
      }
    }
  }
  gsl_matrix* eVec = gsl_matrix_calloc(NDIM4, NDIM4);
  gsl_vector* eVal = gsl_vector_calloc(NDIM4);
  int size = 4 * NDIM4;
  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv(aNew, eVal, eVec, w);
  gsl_eigen_symmv_free(w);
  for (int i = 0; i < rank; i++) {
    eVals->push_back(gsl_vector_get(eVal, i));
  }
  gsl_vector_free(eVal);
  for (int i = 0; i < NDIM4; i++) {
    for (int j = 0; j < NDIM4; j++) {
      // FORTRAN matrices are stored column-major
      eVecs[NDIM4 * i + j] = gsl_matrix_get(eVec, i, j);
    }
  }
  gsl_matrix_free(eVec);
}

inline void Calculate(double spin, double beta2, double beta4, double V0,
                      double R, double A0, double V0S, double A, int Z,
                      int nMax) {
  double SW[2][84] = {};
  double SDW[462] = {};
  int N[NDIM1];
  int L[NDIM1];
  //Matrix elements of Hamiltonian in chosen basis, real & symmetric
  double hamM[NDIM4 * (NDIM4 + 1) / 2];
  double eVecs[NDIM4 * NDIM4];
  double defExpCoef[NDIM4][NDIM1];
  double B[NDIM2];
  double D[NDIM2];
  int LA[NDIM1];
  int IX2[NDIM1];
  int JX2[NDIM1];
  double CNU[NDIM3][NDIM1] = {};
  double E[NDIM3];
  int LK[NDIM3];
  int LKK[NDIM3];
  int KN[7][NDIM3];
  int K1[7];
  double eValsDWS[NDIM4];
  int K2[NDIM4];
  int K3[NDIM4];
  int K4[NDIM4];

  cout << "Entered Calculate from nilsson ns" << endl;

  WoodsSaxon(V0, R, A0, V0S, A, Z, nMax, SW, SDW);

  int II = 0;
  int K = 0;

  int nMin = nMax%2;
  // Diagonalize each L-value separately
  for (int LI = nMin; LI <= nMax; LI += 2) {
    int NIM = II;
    int NI = II + 1;
    int NK = K + 1;

    // Set up quantum numbers of the harmonic oscillator basis
    for (int NN = LI; NN <= nMax; NN += 2) {
      for (int I = 1; I <= 2; I++) {
        N[II] = NN;
        L[II] = LI;
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
    for (int i = NI - 1; i < II; i++) {
      for (int j = NI - 1; j <= i; j++) {
        hamM[KK] = 0.0;
        if (JX2[j] == JX2[i]) {
          int NN = (N[i] / 2) * (N[i] / 2 + 1) * (N[i] / 2 + 2) / 6 +
                   (N[j] / 2) * (N[j] / 2 + 1) / 2 + L[i] / 2;
          hamM[KK] = SW[0][NN] + SW[1][NN] * IX2[j] * (L[j] + LA[j]);
        }
        KK++;
      }
    }
    // inline void Eigen(double* A, int rank, double* eVecs,
    // std::vector<double>* eVals) {
    std::vector<double>* eVals = new std::vector<double>();
    //TODO see if Eigen function returns same values, and in same order when rank != dim
    Eigen(hamM, II - NIM, eVecs, eVals);

    for (int I = NI; I <= II; I++) {
      double eigenValue = (*eVals)[K];
      if (eigenValue <= 10.0) {
        if (K >= NDIM3 - 1) {
          return;
        }
        E[K] = eigenValue;
        LK[K] = L[I - 1];
        K3[K] = JX2[I - 1];
        for (int J = NI; J <= II; J++) {
          int N0 = (II - NIM) * (I - NI) + J - NIM - 1;
          CNU[K][J - 1] = eVecs[N0];
        }
        K++;
      }
    }
    // II = size of harmonic oscillator basis
    // K = size of Woods-Saxon basis for deformed state diagonalization
    if (NK == K + 1) {
      II = NIM;
    }
  }
  if (K == 0) {
    return;
  }

  // TODO Print stuff, line ABOV0161

  int ISX2 = std::abs(2 * spin);
  int ISPIN = (ISX2 + 1) / 2;

  int KK = 0;
  int KKKK = 0;
  int IOM = 1;
  for (int IIOM = 0; IIOM < ISPIN; IIOM++) {
    IOM = -IOM + 2 * (int)(std::pow(-1., IIOM));
    int IIIOM = 4 * (std::abs(IOM) / 4) + 2 - min;
    int KKK = 0;
    for (int I = 0; I < K; I++) {
      LKK[KKK] = LK[I];
      KN[IIOM][KKK] = I;
      double X = 0.0;
      for (int J = 0; J < II; J++) {
        if (L[J] == LKK[KKK]) {
          double cg = utilities::ClebschGordan(2 * L[J], 1, JX2[J],
                                               2 * (LA[J] + IIOM - 1), IX2[J],
                                               2 * (LA[J] + IIOM - 1) + IX2[J]);
          double X1 = CNU[I][J] * cg;
          cg = utilities::ClebschGordan(2 * L[J], 1, JX2[J] - 2 * IX2[J],
                                        2 * (LA[J] + IIOM - 1), IX2[J],
                                        2 * (LA[J] + IIOM - 1) + IX2[J]);
          if (JX2[J] > IIIOM) {
            X1 += CNU[I][J - IX2[J]] * cg;
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
      return;
    }
    K1[IIOM] = KKK;
    if (KKK == 0) {
      break;
    }
    for (int I = 0; I < KKK; I++) {
      for (int J = 0; J < I; J++) {
        B[KK] = 0.0;
        D[KK] = 0.0;
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
                    NI * NJ * (NJ - 1) / 2 + (LI + 1) * NJ + LJ;
              double X = defExpCoef[I][N1] * SDW[NNP - 1] * defExpCoef[J][N2];
              int LP = L[N1];
              int LAP = LA[N1] + IIOM - 1;
              int LL = L[N2];
              int LLA = LA[N2] + IIOM - 1;
              B[KK] +=
                  X * utilities::SphericalHarmonicME(LP, LAP, 2, 0, LL, LLA);
              D[KK] +=
                  X * utilities::SphericalHarmonicME(LP, LAP, 4, 0, LL, LAP);
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
  for (int IIOM = 0; IIOM < ISPIN; IIOM++) {
    IOM = -IOM + 2 * (int)(std::pow(-1., IIOM));
    int KKK = K1[IIOM];
    if (KKK != 0) {
      int NK = 0;
      for (int I = 1; I <= KKK; I++) {
        for (int J = 1; J <= I; J++) {
          hamM[NK] = beta2 * B[KK] + beta4 * D[KK];
          NK++;
          KK++;
        }
        int N0 = KN[IIOM][I];
        hamM[NK] += E[N0];
      }
      std::vector<double>* eVals = new std::vector<double>();
      Eigen(hamM, KKK, eVecs, eVals);
      int N0 = 0;
      for (int MU = 0; MU < KKK; MU++) {
        K2[K] = IOM;
        K3[K] = std::abs(IOM);
        // TODO check
        K4[K] = KKK - MU + 2;
        eValsDWS[K] = hamM[N0];
        for (int J = 0; J < II; J++) {
          defExpCoef[K][J] = 0.0;
          NK = KKK * MU;
          for (int NU = 0; NU < KKK; NU++) {
            int I = KN[IIOM][NU];
            NK++;
            defExpCoef[K][J] += eVecs[NK] * CNU[I][J];
          }
        }
        N0 += MU;
        K++;
      }
    }
  }
}

inline NuclearState CalculateDeformedState(int Z, int A, double beta2,
                                           double beta4) {
  double O;
  double K;
  double parity;

  // TODO

  // special case for 19Ne
  /*struct WFComp s12 = {0.528947, 0, 1, 1};
  struct WFComp d32 = {-0.297834, 2, -1, 1};
  struct WFComp d52 = {0.794676, 2, 1, 1};

  std::vector<WFComp> states = {s12, d32, d52};*/

  // special case for 33Cl
  WFComp d32 = {0.913656, 2, -1, 3};
  WFComp d52 = {-0.406489, 2, 1, 3};

  std::vector<WFComp> states = {d32, d52};

  O = 3. / 2.;
  K = 3. / 2.;
  parity = 1;

  NuclearState nuclearState = {O, K, parity, states};

  return nuclearState;
}
}
#endif
