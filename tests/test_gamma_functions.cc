#include "catch.hpp"

#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_result.h"

#include <cmath>

TEST_CASE("Test real gamma function") {
  REQUIRE(gsl_sf_gamma(1) == 1);
  REQUIRE(gsl_sf_gamma(2) == 1);
  REQUIRE(gsl_sf_gamma(3) == 2);
  REQUIRE(gsl_sf_gamma(0.5) == Approx(1.772).epsilon(0.01));
  REQUIRE(gsl_sf_gamma(1.5) == Approx(0.886).epsilon(0.01));
}

TEST_CASE("Test complex gamma function") {
  gsl_sf_result magn;
  gsl_sf_result phase;
  SECTION("Real value") {
    double zr = 0.5;
    double zi = 0;
    gsl_sf_lngamma_complex_e(zr, zi, &magn, &phase);
    REQUIRE(std::exp(magn.val) == Approx(1.772).epsilon(0.01));
    REQUIRE(phase.val == 0);
  }
  SECTION("Identity for complex value") {
    double zr = 0.5;
    double zi = 10;
    gsl_sf_lngamma_complex_e(zr, zi, &magn, &phase);
    REQUIRE(std::exp(2*magn.val) == Approx(M_PI/std::cosh(M_PI*zi)));
  }
  SECTION("Imaginary value") {
    double zr = 0;
    double zi = 10;
    gsl_sf_lngamma_complex_e(zr, zi, &magn, &phase);
    REQUIRE(std::exp(2*magn.val) == Approx(M_PI/(zi*std::sinh(M_PI*zi))));
  }
}
