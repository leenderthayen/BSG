#ifndef NILSSON_ORBITS
#define NILSSON_ORBITS

#define NDIM1 56
#define NDIM2 820
#define NDIM3 35

namespace nilsson {
  inline void Calculate(double spin, double beta2, double beta4, double V0, double R0, double A0, double V0S, double AMU, int Z, int nMax) {
    double SW[2][84] = {};
    double SDW[462] = {};
    int N[NDIM1];
    int L[NDIM1];
    double AM[NDIM2];
    double R[1600];
    double A[40][NDIM1];
    double B[NDIM2];
    double D[NDIM2];
    double LA[NDIM1];
    double IX2[NDIM1];
    double JX2[NDIM1];
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

    double FACLOG[500];
    int KSELC;

    if (V0 > 0.) {
      Saxon(V0, R0, A0, V0S, AMU, Z, nMax, SW, SDW);
    }
    else if (V0 < 0.) {
      //TODO
    }
    
    int II = 0;
    int K = 0;
    int nMaxP1 = nMax+1;
    int min = nMaxP1/2-nMax/2+1;

    //Diagonalize each L-value separately
    for (int LI=min; LI <= nMaxP1; LI+=2) {
      int NIM = II;
      int NI = II+1;
      int NK = K+1;

      //Set up quantum numbers of the harmonic oscillator basis
      for (int NN=LI; NN <= nMaxP1; NN+=2) {
        for (int I=1; I <= 2; I++) {
          II++;
          N[II] = NN-1;
          L[II] = LI-1;
          LA[II] = 2-I;
          IX2[II] = 2*I-3;
          JX2[II] = 2*L[II]+IX2[II];
          if (L[II] < LA[II]) {
            II--;
          }
        }
      }
      if (II == NIM) {
        break;
      }

      //Set up matrix
      int KK = 0;
      for (int I = NI; I <= II; I++) {
        for (int J = NI; J <= I; J++) {
          KK++;
          AM[KK] = 0.0;
          if (!(JX2[J] != JX2[I])) {
            NN = (N[I]/2)*(N[I]/2+1)*(N[I]/2+2)/6+(N[J]/2)*(N[J]/2+1)/2+L[I]/2+1;
            AM[KK] = SW[1][NN]+SW[2][NN]*IX2[J]*(L[J]+LA[J]);
          }
        }
      }
      Eigen(AM, R, II-NIM, 0);

      for (int I = NI; I <= II; I++) {
        int N0 = (I-NIM)*(I-NIM+1)/2;
        if (AM[N0] <= 10.0) {
          K++;
          if (K >= NDIM3) {
            return;
          }
          E[K] = AM[N0];
          LK[K] = L[I];
          K3[K] = KX2[I];
          for (int J = NI; J <= II; J++) {
            N0 = (II-NIM)*(I-NI)+J-NIM;
            CNU[K][J] = R[N0];
          }
        }
      }
      //II = size of harmonic oscillator basis
      //K = size of Woods-Saxon basis for deformed state diagonalization
      if (NK == K+1) {
        II = NIM;
      }
    }
    if (K == 0) {
      return;
    }

    //TODO Print stuff, line ABOV0161

    //Initialize Clebsch-Gordan Faclog table
    KSELC = 1;
    FACLOG[1] = 0.0;
    FACLOG[2] = 0.0;
    double X  = 0.0;
    for (int i=3; i <= 50; i++) {
      X++;
      FACLOG[i] = FACLOG[i-1] + std::log(X);
    }
    int ISX2 = std::abs(2*spin);
    int ISPIN = (ISX2+1)/2;
    
  }
}
#endif
