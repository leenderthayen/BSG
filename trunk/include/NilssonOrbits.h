#ifndef NILSSON_ORBITS
#define NILSSON_ORBITS

#include "Utilities.h"
#include "Constants.h"

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#define NDIM1 56
#define NDIM2 820
#define NDIM3 35
#define NDIM4 40

namespace nilsson {

//Here dO stands for double Omega
//s is +-1 depending of whether it is j=l+-1/2
struct WFComp {double C; int l, s, dO; };

struct NuclearState { double O, K; int parity; std::vector<WFComp> states; };

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
  double vNorm = 2./std::pow(M_PI, 0.25)*eMinus;
  int k = (n+1)/2-n/2+1;
  if (k > 1) {
    vNorm *= std::sqrt(2./3.);
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

inline void WoodsSaxon(double V0, double R, double A0, double V0S, int A, int Z, int nMax, double SW[2][84], double SDW[462]) {
  double FINT[3] = {};
  double S[3][4] = {};
  int IP = 2;
  int IW = 3;
  int II = 0;
  int JJ = 0;
  int NDX = 300;
  int N = NDX/4-2;
  double A = std::pow(A, 1./3.);
  double AO = A0;
  //Compton wavelength of the pion=  1.4133 fm
  double pionComp = hbar/pionMasskeV/speedOfLight;
  double VOS = V0S*pionComp*pionComp;
  double TWONU = 0.154666*std::sqrt(V0)/R;
  double T = TWONU*20.899*A/(A-1.);
  //e^2 = 1.439978 MeV*fm
  double ZZ = Z*1.439978;
  double ZCONST = 1.5*ZZ/R;
  double ZIN = 0.5*ZZ/std::pow(R, 3.);
  double DX = (R+8.*AO)/NDX;
  double DFO = std::exp(DX/AO);
  int NMAXP1 = nMax+1;
  //equal to 1 when even, 2 when odd
  int NMIN = NMAXP1/2-(NMAXP1-1)/2+1;
  for (int NNI = NMIN; NNI <= NMAXP1; NNI+=2) {
    for (int NNJ = NMIN; NNJ <= NNI; NNJ+=2) {
      for (int LLI = NMIN; LLI <= NNI; LLI+=2) {
        for (int LLJ = NMIN; LLJ <= NNJ; LLJ+=2) {
          int NI = NNI-1;
          int NJ = NNJ-1;
          int LI = LLI-1;
          int LJ = LLJ-1;

          int KK = 0;
          if (LI != LJ) {
            KK = 2;
          }
          double r = DX;
          double FO = DFO/std::exp(R/AO);
          int M = 0;
          while(M/N >= 0) {
            if (M == 0) {
              M = -N-2;
            }
            for (int I = 0; I < 4; I++) {
              double R2 = r*r*TWONU;
              double PSI2 = V(NI, LI, R2)*V(NJ, LJ, R2)*std::pow(R2, (LI+LJ)/2)*R*R/std::exp(R2);
              double F = 1./(1.+FO);
              double DFDR = -FO*F*F/AO;
              if (r <= R) {
                S[0][I] = -PSI2*(V0*F+T*R2+ZIN*r*r-ZCONST);
              } else {
                S[0][I] = -PSI2*(V0*F+T*R2-ZZ/r);
              }
              S[1][I] = PSI2*DFDR/r;
              S[2][I] = PSI2*DFDR;
              FO *= DFO;
              r += DX;
              for (int j = KK; j < 3; j++) {
                FINT[j] += S[j][I];
              }
            }
            M += 2;
            if (M/N > 0) {
              r += 4.*DX;
              DX = -DX;
              FO *= std::pow(DFO, 4.);
              DFO = 1./DFO;
            }
          }
          for (int i = 1; i < 4; i++) {
            for (int l = i; l < 4; l++) {
              //TODO look at K
              int K = 3+i-l;
              for (int j = KK; j < 3; j++) {
                S[j][K] -= S[j][K-1];
              }
            }
          }
          for (int j = KK; j < 3; j++) {
            FINT[j] += -S[j][0]/2.+S[j][1]/12.-S[j][2]/24.+19.*S[j][3]/720.;
          }
          //TODO
          if (DX <= 0) {
            DX = -DX;
            DFO = 1./DFO;
          }
          double X = VNORM(NI, LI)*VNORM(NJ, LJ);
          if (KK != 2) {
            SW[0][II] = FINT[0]*X*DX*std::pow(TWONU, 1.5);
            if (NI == NJ) {
              SW[0][II] += 2.*T*(NI+1.5);
            }
            SW[1][II] = FINT[1]*X*DX*std::pow(TWONU, 1.5)*VOS;
            II++;
          }
          SDW[JJ] = FINT[2]*X*DX*std::pow(TWONU, 1.5)*V0*RO;
          JJ++;
        }
      }
    }
  }
}

//Calculate eigenvalues & vectors for a real symmetric FORTRAN matrix A
inline void Eigen(double* A, int rank, double* eVecs, std::vector<double>* eVals) {
  gsl_matrix * aNew = gsl_matrix_alloc(NDIM4, NDIM4);
  for (int i = 0; i < NDIM4; i++) {
    for (int j = i; j < NDIM4; j++) {
      //FORTRAN matrices are stored column-major
      gsl_matrix_set(aNew, i, j, A[NDIM*i+j]);
      if (i != j) {
        gsl_matrix_set(aNew, j, i, A[NDIM4*i+j]);
      }
    }
  }
  gsl_matrix* eVec = gsl_matrix_calloc(NDIM4, NDIM4);
  gsl_vector* eVal = gsl_vector_calloc(NDIM4);
  int size = 4*NDIM4;
  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv(aNew, eVal, eVec, w);
  gsl_eigen_symmv_free(w);
  for (int i = 0; i < rank; i++) {
    eVals->push_back(gsl_vector_get(eVal, i));
  }
  gsl_vector_free(eVal);
  for (int i = 0; i < NDIM4; i++) {
    for (int j = 0; j < NDIM4; j++) {
      //FORTRAN matrices are stored column-major
      eVecs[NDIM4*i+j] = gsl_matrix_get(eVec, i, j);
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
  double AM[NDIM4*(NDIM4+1)/2];
  double R[NDIM4*NDIM4];
  double A[NDIM4][NDIM1];
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
  double ED[NDIM4];
  int K2[NDIM4];
  int K3[NDIM4];
  int K4[NDIM4];

  WoodsSaxon(V0, R0, A0, V0S, A, Z, nMax, SW, SDW);

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
        N[II] = NN - 1;
        L[II] = LI - 1;
        LA[II] = 2 - I;
        IX2[II] = 2 * I - 3;
        JX2[II] = 2 * L[II] + IX2[II];
        if (L[II] < LA[II]) {
          II--;
        }
        II++;
      }
    }
    if (II == NIM) {
      break;
    }

    // Set up matrix
    int KK = 0;
    for (int I = NI-1; I < II; I++) {
      for (int J = NI-1; J < I; J++) {
        AM[KK] = 0.0;
        if (!(JX2[J] != JX2[I])) {
          int NN = (N[I] / 2) * (N[I] / 2 + 1) * (N[I] / 2 + 2) / 6 +
               (N[J] / 2) * (N[J] / 2 + 1) / 2 + L[I] / 2;
          AM[KK] = SW[0][NN] + SW[1][NN] * IX2[J] * (L[J] + LA[J]);
        }
        KK++;
      }
    }
    //inline void Eigen(double* A, int rank, double* eVecs, std::vector<double>* eVals) {
    std::vector<double> * eVals = new std::vector<double>();
    Eigen(AM, II - NIM, R, eVals);

    for (int I = NI; I <= II; I++) {
      double eigenValue = (*eVals)[K];
      if (eigenValue <= 10.0) {
        if (K >= NDIM3-1) {
          return;
        }
        E[K] = eigenValue;
        LK[K] = L[I-1];
        K3[K] = JX2[I-1];
        for (int J = NI; J <= II; J++) {
          N0 = (II - NIM) * (I - NI) + J - NIM - 1;
          CNU[K][J-1] = R[N0];
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
          A[KKK][J] = X1;
          X += X1 * X1;
        }
      }
      if (X < 0.1) {
        KKK--;
      }
      KKK++;
    }
    KKKK+=KKK;
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
              double X = A[I][N1]*SDW[NNP-1]*A[J][N2];
              int LP = L[N1];
              int LAP = LA[N1]+IIOM-1;
              int LL = L[N2];
              int LLA = LA[N2]+IIOM-1;
              B[KK]+=X*utilities::SphericalHarmonicME(LP, LAP, 2, 0, LL, LLA);
              D[KK]+=X*utilities::SphericalHarmonicME(LP, LAP, 4, 0, LL, LAP);
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
    IOM = -IOM+2*(int)(std::pow(-1., IIOM));
    int KKK = K1[IIOM];
    if (KKK != 0) {
      int NK = 0;
      for (int I = 1; I <= KKK; I++) {
        for (int J = 1; J <= I; J++) {
          AM[NK] = beta2*B[KK]+beta4*D[KK];
          NK++;
          KK++;
        }
        int N0 = KN[IIOM][I];
        AM[NK]+=E[N0];
      }
      Eigen(AM, R, KKK, 0);
      int N0 = 0;
      for (int MU = 0; MU < KKK; MU++) {
        K2[K] = IOM;
        K3[K] = std::abs(IOM);
        //TODO check
        K4[K] = KKK-MU+2;
        ED[K] = AM[N0];
        for (int J = 0; J < II; J++) {
          A[K][J] = 0.0;
          NK = KKK*MU;
          for (int NU = 0; NU < KKK; NU++) {
            int I = KN[IIOM][NU];
            NK++;
            A[K][J]+R[NK]*CNU[I][J];
          }
        }
        N0+=MU;
        K++;
      }
    }
  }
}

inline NuclearState CalculateDeformedState(int Z, int A, double beta2, double beta4) {
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
  WFComp d32 = {0.913656, 2, -1, 3};
  WFComp d52 = {-0.406489, 2, 1, 3};

  std::vector<WFComp> states = {d32, d52};

  O = 3./2.;
  K = 3./2.;
  parity = 1;

  NuclearState nuclearState = {O, K, parity, states};

  return nuclearState;
}
}
#endif
