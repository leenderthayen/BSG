#ifndef NILSSON_ORBITS
#define NILSSON_ORBITS

#include "utilities.h"
#include "gsl/gsl_eigen.h"

#define NDIM1 56
#define NDIM2 820
#define NDIM3 35

namespace nilsson {

inline double V(int n, int l, double x) {
  int m = (n-l)/2;
  double V = 1.0;
  double F = 1.0;
  int mK = 1;
  int mM = 1;
  if (m > 0) {
    for (int k = 1; k <= m; k++) {
      mK *= m+1-k;
      F *= -2.*x/k;
      mM *= 2*n-4*m+1+2*k;
      V += F*mK/mM;
    }
    return V;
  }
  return 0;
}

inline double VNORM(int n, int l) {
  int minus = (n-l)/2;
  int eMinus = 1+4*(minus/2)-2*minus;
  double vNorm = 1.50225*eMinus;
  int k = (n+1)/2-n/2+1;
  if (k > 1) {
    vNorm *= 0.816496;
  }
  if (n > 1) {
    for (int j = k; j <= n; j+=2) {
      vNorm *= std::sqrt((double)(j+k+1.)/(j-k+2.));
    }
  }
  if (l > 1) {
    for (int j = k; j <= l; j+=2) {
      vNorm *= 2.*std::sqrt((double)((j+n+2)*(n-j+1))/((2.*j+3.)*(2.*j+1.)));
    }
  }
  return vNorm;
}

inline void WoodsSaxon(double V0, double R0, double A0, double V0S, double AMU, int Z, int nMax, double SW[2][84], double SDW[462]) {
  double FINT[3] = {};
  double S[3][4] = {};
  int IP = 2;
  int IW = 3;
  int II = 0;
  int JJ = 0;
  int NDX = 300;
  int N = NDX/4-2;
  double A = std::pow(AMU, 1./3.);
  double RO = A*std::abs(R0);
  double AO = A0;
  double pionComp = 1.4133;
  double VOS = V0S*pionComp*pionComp;
  double TWONU = 0.154666*std::sqrt(V0)/RO;
  double T = TWONU*20.899*AMU/(AMU-1.);
  double ZZ = Z*1.439978;
  double ZCONST = 1.5*ZZ/RO;
  double ZIN = 0.5*ZZ/std::pow(RO, 3.);
  double DX = (RO+8.*AO)/NDX;
  double DFO = std::exp(DX/AO);
  int NMAXP1 = nMax+1;
  int NMIN = NMAXP1/2-(NMAXP1-1)/2+1;
  for (int NNI = NMIN; NNI <= NMAXP1; NNI+=2) {
    for (int NNJ = NMIN; NNJ <= NNI; NNJ+=2) {
      for (int LLI = NMIN; LLI <= NNI; LLI+=2) {
        for (int LLJ = NMIN; LLJ <= NNJ; LLJ+=2) {
          int NI = NNI-1;
          int NJ = NNJ-1;
          int LI = LLI-1;
          int LJ = LLJ-1;

          int KK = 1;
          if (LI != LJ) {
            KK = 3;
          }
          double R = DX;
          double FO = DFO/std::exp(RO/AO);
          for (int i = 0; i < 3; i++) {
            FINT[i] = 0.0;
          }
          int M = 0;
          while(M/N >= 0) {
            if (M == 0) {
              M = -N-2;
            }
            for (int I = 1; I <= 4; I++) {
              double R2 = R*R*TWONU;
              double PSI2 = V(NI, LI, R2)*V(NJ, LJ, R2)*std::pow(R2, (LI+LJ)/2)*R*R/std::exp(R2);
              double F = 1./(1.+FO);
              double DFDR = -FO*F*F/AO;
              if (R <= RO) {
                S[1][I] = -PSI2*(V0*F+T*R2+ZIN*R*R-ZCONST);
              } else {
                S[1][I] = -PSI2*(V0*F+T*R2-ZZ/R);
              }
              S[2][I] = PSI2*DFDR/R;
              S[3][I] = PSI2*DFDR;
              FO *= DFO;
              R += DX;
              for (int j = KK; j <= 3; j++) {
                FINT[j] += S[j][I];
              }
            }
            M += 2;
            if (M/N > 0) {
              R += 4.*DX;
              DX = -DX;
              FO *= std::pow(DFO, 4.);
              DFO = 1./DFO;
            }
          }
          for (int i = 2; i <= 4; i++) {
            for (int l = i; l <= 4; l++) {
              int K = 4+i-l;
              for (int j = KK; j <= 3; j++) {
                S[j][K] -= S[j][K-1];
              }
            }
          }
          for (int j = KK; j <= 3; j++) {
            FINT[j] += -S[j][1]/2.+S[j][2]/12.-S[j][3]/24.+19.*S[j][4]/720.;
          }
          //TODO
          if (DX <= 0) {
            DX = -DX;
            DFO = 1./DFO;
          }
          double X = VNORM(NI, LI)*VNORM(NJ, LJ);
          if (KK != 3) {
            II++;
            SW[1][II] = FINT[1]*X*DX*std::pow(TWONU, 1.5);
            if (NI == NJ) {
              SW[1][II] += 2.*T*(NI+1.5);
            }
            SW[2][II] = FINT[2]*X*DX*std::pow(TWONU, 1.5)*VOS;
          }
          JJ++;
          SDW[JJ] = FINT[3]*X*DX*std::pow(TWONU, 1.5)*V0*RO;
        }
      }
    }
  }
}

inline void Eigen(double A[], double R[], int N, int MV) {

}

inline void Calculate(double spin, double beta2, double beta4, double V0,
                      double R0, double A0, double V0S, double AMU, int Z,
                      int nMax) {
  double SW[2][84] = {};
  double SDW[462] = {};
  int N[NDIM1];
  int L[NDIM1];
  double AM[NDIM2];
  double R[1600];
  double A[40][NDIM1];
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
  double ROTO[7];
  double RROTO[7];
  double ROTI[7];
  double ED[40];
  int K2[40];
  int K3[40];
  int K4[40];

  if (V0 > 0.) {
    WoodsSaxon(V0, R0, A0, V0S, AMU, Z, nMax, SW, SDW);
  } else if (V0 < 0.) {
    // TODO
  }

  int II = 0;
  int K = 0;
  int nMaxP1 = nMax + 1;
  int min = nMaxP1 / 2 - nMax / 2 + 1;

  // Diagonalize each L-value separately
  for (int LI = min; LI <= nMaxP1; LI += 2) {
    int NIM = II;
    int NI = II + 1;
    int NK = K + 1;

    // Set up quantum numbers of the harmonic oscillator basis
    for (int NN = LI; NN <= nMaxP1; NN += 2) {
      for (int I = 1; I <= 2; I++) {
        II++;
        N[II] = NN - 1;
        L[II] = LI - 1;
        LA[II] = 2 - I;
        IX2[II] = 2 * I - 3;
        JX2[II] = 2 * L[II] + IX2[II];
        if (L[II] < LA[II]) {
          II--;
        }
      }
    }
    if (II == NIM) {
      break;
    }

    // Set up matrix
    int KK = 0;
    for (int I = NI; I <= II; I++) {
      for (int J = NI; J <= I; J++) {
        KK++;
        AM[KK] = 0.0;
        if (!(JX2[J] != JX2[I])) {
          int NN = (N[I] / 2) * (N[I] / 2 + 1) * (N[I] / 2 + 2) / 6 +
               (N[J] / 2) * (N[J] / 2 + 1) / 2 + L[I] / 2 + 1;
          AM[KK] = SW[1][NN] + SW[2][NN] * IX2[J] * (L[J] + LA[J]);
        }
      }
    }
    Eigen(AM, R, II - NIM, 0);

    for (int I = NI; I <= II; I++) {
      int N0 = (I - NIM) * (I - NIM + 1) / 2;
      if (AM[N0] <= 10.0) {
        K++;
        if (K >= NDIM3) {
          return;
        }
        E[K] = AM[N0];
        LK[K] = L[I];
        K3[K] = JX2[I];
        for (int J = NI; J <= II; J++) {
          N0 = (II - NIM) * (I - NI) + J - NIM;
          CNU[K][J] = R[N0];
        }
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
  for (int IIOM = 1; IIOM <= ISPIN; IIOM++) {
    IOM = -IOM - 2 * (int)(std::pow(-1., IIOM));
    int IIIOM = 4 * (std::abs(IOM) / 4) + 2 - min;
    int KKK = 0;
    for (int I = 1; I <= K; I++) {
      LKK[KKK] = LK[I];
      KN[IIOM][KKK] = I;
      double X = 0.0;
      for (int J = 1; J <= II; J++) {
        if (L[J] == LKK[KKK]) {
          double cg = utilities::ClebschGordan(2 * L[J], 1, JX2[J],
                                               2 * (LA[J] + IIOM - 1), IX2[J],
                                               2 * (LA[J] + IIOM - 1) + IX2[J]);
          // TODO Check if coefficient I is correct
          double X1 = CNU[I][J] * cg;
          cg = utilities::ClebschGordan(2 * L[J], 1, JX2[J] - 2 * IX2[J],
                                        2 * (LA[J] + IIOM - 1), IX2[J],
                                        2 * (LA[J] + IIOM - 1) + IX2[J]);
          if (JX2[J] > IIIOM) {
            X1 += CNU[I][J - IX2[J]] * cg;
          }
          A[KKK][J] = X1;
          X += X1 * X1;
        }
      }
      if (X < 0.1) {
        KKK--;
      }
    }
    KKKK+=KKK;
    if (KKKK > 40) {
      return;
    }
    K1[IIOM] = KKK;
    if (KKK == 0) {
      break;
    }
    for (int I = 1; I <= KKK; I++) {
      for (int J = 1; J <= I; J++) {
        KK++;
        B[KK] = 0.0;
        D[KK] = 0.0;
        for (int N1 = 1; N1 <= II; N1++) {
          for (int N2 = 1; N2 <= II; N2++) {
            if (!(L[N1] != LKK[I] || L[N2] != LKK[J] || (L[N1]-L[N2])/5 != 0 || LA[N1] != LA[N2])) {
              int NI = N[N1]/2+1;
              int NJ = N[N2]/2+1;
              int LI = L[N1]/2+1;
              int LJ = L[N2]/2+1;
              if (NI < NJ) {
                int NN = NI;
                NI = NJ;
                NJ = NN;
                NN = LI;
                LI = LJ;
                LJ = NN;
              }
              int NNP = NI*(NI-1);
              NNP = (3*NNP*NNP+2*NNP*(2*NI-1))/24+NI*NJ*(NJ-1)/2+(LI+1)*NJ+LJ;
              double X = A[I][N1]*SDW[NNP]*A[J][N2];
              int LP = L[N1];
              int LAP = LA[N1]+IIOM-1;
              int LL = L[N2];
              int LLA = LA[N2]+IIOM-1;
              B[KK]+=X*utilities::SphericalHarmonicME(LP, LAP, 2, 0, LL, LLA);
              D[KK]+=X*utilities::SphericalHarmonicME(LP, LAP, 4, 0, LL, LAP);
            }
          }
        }
      }
    }
  }

  K = 0;
  KK = 0;
  IOM = 1;
  for (int IIOM = 1; IIOM <= ISPIN; IIOM++) {
    IOM = -IOM-2*(int)(std::pow(-1., IIOM));
    int KKK = K1[IIOM];
    if (KKK != 0) {
      int NK = 0;
      for (int I = 1; I <= KKK; I++) {
        for (int J = 1; J <= I; J++) {
          NK++;
          KK++;
          AM[NK] = beta2*B[KK]+beta4*D[KK];
        }
        int N0 = KN[IIOM][I];
        AM[NK]+=E[N0];
      }
      Eigen(AM, R, KKK, 0);
      int N0 = 0;
      for (int MU = 1; MU <= KKK; MU++) {
        N0+=MU;
        K++;
        K2[K] = IOM;
        K3[K] = std::abs(IOM);
        K4[K] = KKK-MU+1;
        ED[K] = AM[N0];
        for (int J = 1; J <= II; J++) {
          A[K][J] = 0.0;
          NK = KKK*(MU-1);
          for (int NU = 1; NU <= KKK; NU++) {
            int I = KN[IIOM][NU];
            NK++;
            A[K][J]+R[NK]*CNU[I][J];
          }
        }
      }
    }
  }
}

inline utilities::NuclearState CalculateDeformedState(int Z, int A, double beta2, double beta4) {
  double O;
  double K;
  double parity;

  //TODO

  //special case for 19Ne
  /*struct WFComp s12 = {0.528947, 0, 1, 1};
  struct WFComp d32 = {-0.297834, 2, -1, 1};
  struct WFComp d52 = {0.794676, 2, 1, 1};
  
  std::vector<WFComp> states = {s12, d32, d52};*/

  //special case for 33Cl
  utilities::WFComp d32 = {0.913656, 2, -1, 3};
  utilities::WFComp d52 = {-0.406489, 2, 1, 3};

  std::vector<utilities::WFComp> states = {d32, d52};

  O = 3./2.;
  K = 3./2.;
  parity = 1;

  utilities::NuclearState nuclearState = {O, K, parity, states};

  return nuclearState;
}
}
#endif
