#include "catch.hpp"

#include "BSG/SpectralFunctions.h"

TEST_CASE("Spectral functions: phase space") {
  double W0 = 10;
  REQUIRE(BSG::SpectralFunctions::PhaseSpace(1, W0) == 0);
  REQUIRE(BSG::SpectralFunctions::PhaseSpace(W0, W0) == 0);
  REQUIRE(BSG::SpectralFunctions::PhaseSpace(5, W0) == Approx(612.37));
}

TEST_CASE("Spectral functions: Fermi function") {
  double R = 0.01;
  REQUIRE(BSG::SpectralFunctions::FermiFunction(1, 0, R, 1) == 1);
  REQUIRE(BSG::SpectralFunctions::FermiFunction(10, 0, R, 1) == 1);

  //Pretty strong dependence on R, so don't aim higher than this to make sure it works
  double eps = 5e-3;

  //Test cases come from The Landolt-Boernstein tables.
  //Note that there is a subtle difference here, as these values include the L0 correction
  SECTION("Low Z") {
    int Z = 1;
    R = 0.005;
    REQUIRE(BSG::SpectralFunctions::FermiFunction(1.005, Z, R, 1) == Approx(1.2487).epsilon(eps));
    REQUIRE(BSG::SpectralFunctions::FermiFunction(10, Z, R, 1) == Approx(1.0233).epsilon(eps));
  }
  SECTION("High Z") {
    int Z = 50;
    R = 0.015;
    REQUIRE(BSG::SpectralFunctions::FermiFunction(1.00499, Z, R, 1) == Approx(55.613/1.0153).epsilon(eps));
    REQUIRE(BSG::SpectralFunctions::FermiFunction(10.05, Z, R, 1) == Approx(3.9468/0.9572).epsilon(eps));
    REQUIRE(BSG::SpectralFunctions::FermiFunction(50.01, Z, R, 1) == Approx(2.4679/0.7498).epsilon(eps));
  }
}
