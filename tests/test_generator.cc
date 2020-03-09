#include "catch.hpp"

#include "BSG/Generator.h"
#include "NHL/Containers.h"

TEST_CASE("Basic construction") {
  std::string iniFile = "/Users/leenderthayen/git/BSG/tests/45Ca.ini";
  BSG::Generator gen = BSG::Generator();

  bool success = gen.InitialiseTransitionFromFile(iniFile);

  REQUIRE(success == true);

  BSG::BetaParams bp = gen.GetBetaParams();

  REQUIRE(bp.Zi == 20);
  REQUIRE(bp.A == 45);
  REQUIRE(bp.Zf == 21);
  REQUIRE(bp.betaType == NHL::BETA_MINUS);
  REQUIRE(bp.decayType == NHL::BetaDecayType::GAMOWTELLER);
  REQUIRE(bp.R == Approx(0.011).epsilon(0.1));
  REQUIRE(bp.mixingRatio == 0.);
  REQUIRE(bp.aNeg[0] != 0.);
  REQUIRE(bp.W0 == Approx(1.5).epsilon(0.01));

  std::vector< std::vector<double> >* spectrum = gen.CalculateSpectrum();

  REQUIRE((*spectrum)[0][0] != 0);
  REQUIRE((*spectrum)[0][1] != 0);
  REQUIRE((*spectrum)[1][0] != 0);
}
