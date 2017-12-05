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

/**
 * Shell model occupations in the jj coupling scheme
 * each element consists of 4 numbers.
 * Example: 1d5/2 state becomes {1, 2, 1, 14}
 * Here the second number of the orbital angular momentum
 * the third is +1 for the l+1/2 component and -1 otherwise
 * the last is the total occupation when all shells up to that one
 * is filled
 */
static int smOccupation[104] = {
    1, 0, 1,  2,   1, 1, 1,  6,   1, 1, -1, 8,   1, 2, 1,  14,  2, 0, 1,  16,
    1, 2, -1, 20,  1, 3, 1,  28,  2, 1, 1,  32,  1, 3, -1, 38,  2, 1, -1, 40,
    1, 4, 1,  50,  1, 4, -1, 58,  2, 2, 1,  64,  2, 2, -1, 68,  3, 0, 1,  70,
    1, 5, 1,  82,  1, 5, -1, 92,  2, 3, 1,  100, 2, 3, -1, 106, 3, 1, 1,  110,
    3, 1, -1, 112, 1, 6, 1,  126, 2, 4, 1,  136, 3, 2, 1,  142, 1, 6, -1, 154};

/**
 * Spectroscopic names corresponding to orbital angular momentum
 */
static char spectroNames[7] = {'s', 'p', 'd', 'f', 'g', 'h', 'i'};

/**
 * Calculate the factorial n!
 * 
 * @param n integer
 */
inline int Factorial(int n) {
  if (n <= 0) {
    return 1;
  } else {
    return n * Factorial(n - 1);
  }
}

/**
 * Double factorial function n!!
 * 
 * @param n integer
 */
inline int DoubleFactorial(int n) {
  if (n <= 0) {
    return 1;
  } else {
    return n * DoubleFactorial(n - 2);
  }
}

/**
 * Calculate the Clebsch-Gordan coefficient
 * 
 * @param two_ja double of first spin
 * @param two_jb double of second spin
 * @param two_jc double of resultant spin
 * @param two_ma double of z projection of first spin
 * @param two_mb double of z projection of second spin
 * @param two_mc double of z projection of resultant spin
 */
inline double ClebschGordan(int two_ja, int two_jb, int two_jc, int two_ma,
                            int two_mb, int two_mc) {
  double result = 0.0;
  if (two_ja >= 0 && two_jb >= 0 && two_jc >= 0) {
    result =
        std::pow(-1., (two_ja - two_jb + two_mc) / 2) *
        std::sqrt(two_jc + 1.0) *
        gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, -two_mc);
  }
  return result;
}

// Function to calculate the matrix elements of spherical harmonics
/*                       M
           < L', LA' | Y  | L, LA >
                         L
*/
/**
 * Calculate the matrix element of a spherical harmonic
 * @f[\langle J_f M_f | Y^M_L | J_i M_i \rangle @f]
 * 
 * @param lP final spin
 * @param laP z-projection of the final spin
 * @param l rank of spherical tensor of the operator
 * @param m z-projection of the rank of the operator
 * @param ll initial spin
 * @param la z-projection of the initial spin
 */
inline double SphericalHarmonicME(int lP, int laP, int l, int m, int ll,
                                  int la) {
  double result = 0.;
  if (laP == m + la) {
    double x = (2. * l + 1.) * (2. * ll + 1.) / (12.566375 * (2. * lP + 1.));
    double cg =
        ClebschGordan(2 * ll, 2 * l, 2 * lP, 2 * la, 2 * m, 2 * (la + m));
    result = std::sqrt(x) * cg;
    cg = ClebschGordan(2 * ll, 2 * l, 2 * lP, 0, 0, 0);
    result *= cg;
  }
  return result;
}

/**
 * Calculate the occupation numbers of the shell model orbitals assuming the jj-coupling scheme
 * 
 * @param N number of nucleons
 */
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

/**
 * Class for Lagrange interpolation between a set of 3 (x,y) points
 */
class Lagrange {
 public:
  Lagrange(double*, double*);
  ~Lagrange(){};
  double GetValue(double);

 private:
  double xC[3];
  double yC[3];
};

/**
 * Perform Simpson integration
 * 
 * @param x array of x values
 * @param y array of y values
 * @param size size of the arrays
 */
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

inline double Simpson(std::vector<std::vector<double> > values) {
  double result = 0.;
  for (int i = 0; i < values.size()-2; i++) {
    double xN[] = {values[i][0], values[i][0], values[i][0]};
    double yN[] = {values[i][1], values[i][1], values[i][1]};
    double h = (xN[2] - xN[0]) / 2.;
    Lagrange* l = new Lagrange(xN, yN);
    result += 1. / 3. * h * (values[i][1] + 4. * l->GetValue(xN[0] + h) + values[i + 2][1]);
    if (result != result) result = 0.;
  }
  return result;
}

/**
 * Perform trapezoid integration
 * 
 * @param x array of x values
 * @param y array of y values
 * @param array size
 */
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
