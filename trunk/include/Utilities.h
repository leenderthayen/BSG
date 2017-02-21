#ifndef UTILITIES
#define UTILITIES

// standard classes
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>

#include "Constants.h"

#include "gsl/gsl_sf_coupling.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multifit_nlin.h"

namespace utilities {

using std::cout;
using std::endl;

static int smOccupation[104] = {1, 0, 1, 2, 1, 1, 1, 6,  1, 1, -1, 8, 1, 2, 1, 14,  2, 0, 1, 16,
                        1, 2, -1, 20, 1, 3, 1, 28, 2, 1, 1, 32, 1, 3, -1, 38, 2, 1, -1, 40,
                        1, 4, 1, 50, 1, 4, -1, 58, 2, 2, 1, 64, 2, 2, -1, 68, 3, 0, 1, 70,
                        1, 5, 1, 82,  1, 5, -1, 92,  2, 3, 1, 100, 2, 3, -1, 106, 3, 1, 1, 110,
                        3, 1, -1, 112, 1, 6, 1, 126, 2, 4, 1, 136, 3, 2, 1, 142, 1, 6, -1, 154};

inline int Factorial(int n) {
  if (n <= 0) {
    return 1;
  } else {
    return n * Factorial(n - 1);
  }
}

inline int DoubleFactorial(int n) {
  if (n <= 0) {
    return 1;
  } else {
    return n * DoubleFactorial(n - 2);
  }
}

inline double ClebschGordan(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc) {
  return std::pow(-1., (two_ja-two_jb+two_mc)/2)*std::sqrt(two_jc+1)*gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, -two_mc);
}

// Function to calculate the matrix elements of spherical harmonics
/*                       M
           < L', LA' | Y  | L, LA >
                         L

*/
inline double SphericalHarmonicME(int lP, int laP, int l, int m, int ll, int la) {
  double result = 0.;
  if (laP == m+la) {
    double x = (2.*l+1.)*(2.*ll+1.)/(12.566375*(2.*lP+1.));
    double cg = ClebschGordan(2*ll, 2*l, 2*lP, 2*la, 2*m, 2*(la+m));
    result = x*cg;
    cg = ClebschGordan(2*ll, 2*l, 2*lP, 0, 0, 0);
    result*=cg;
  }
  return result;
}

inline std::vector<int> GetOccupationNumbers(int N) {
  std::vector<int> occNumbers;
  occNumbers.push_back(1);
  occNumbers.push_back(0);
  occNumbers.push_back(1);
  occNumbers.push_back(std::min(N, 2));
  for (int i = 4; i < sizeof(smOccupation) / sizeof(*smOccupation); i += 4) {
    if (N > smOccupation[i - 1]) {
      occNumbers.push_back(smOccupation[i]);
      occNumbers.push_back(smOccupation[i + 1]);
      occNumbers.push_back(smOccupation[i + 2]);
      occNumbers.push_back(std::min(N - smOccupation[i - 1],
                               smOccupation[i + 3] - smOccupation[i - 1]));
    }
  }
  return occNumbers;
}


class Lagrange {
 public:
  Lagrange(double*, double*);
  ~Lagrange(){};
  double GetValue(double);

 private:
  double xC[3];
  double yC[3];
};

inline double Simpson(double x[], double y[], int size) {
  double result = 0.;
  for (int i = 0; i < size - 2; i += 2) {
    double xN[] = {x[i], x[i + 1], x[i + 2]};
    double yN[] = {y[i], y[i + 1], y[i + 2]};
    double h = (xN[2] - xN[0]) / 2.;
    Lagrange* l = new Lagrange(xN, yN);
    result += 1. / 3. * h * (y[i] + 4. * l->GetValue(xN[0] + h) + y[i + 2]);
    if (result != result) result = 0.;
  }
  // TODO add result for last values
  return result;
}

inline double Trapezoid(double x[], double y[], int size) {
  double result = 0.;
  for (int i = 0; i < size - 1; i++) {
    double h = x[i + 1] - x[i];
    result += h / 2. * (y[i + 1] + y[i]);
  }
  return result;
}
}
#endif
