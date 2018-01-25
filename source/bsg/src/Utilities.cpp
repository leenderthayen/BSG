#include "Utilities.h"

utilities::Lagrange::Lagrange(double* x, double* y) {
  xC[0] = x[0];
  xC[1] = x[1];
  xC[2] = x[2];
  yC[0] = y[0];
  yC[1] = y[1];
  yC[2] = y[2];
}

double utilities::Lagrange::GetValue(double x) {
  double first =
      yC[0] * (x - xC[1]) * (x - xC[2]) / (xC[0] - xC[1]) / (xC[0] - xC[2]);
  double second =
      yC[1] * (x - xC[0]) * (x - xC[2]) / (xC[1] - xC[0]) / (xC[1] - xC[2]);
  double third =
      yC[2] * (x - xC[0]) * (x - xC[1]) / (xC[2] - xC[0]) / (xC[2] - xC[1]);

  return first + second + third;
}
