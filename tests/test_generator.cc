#include "catch.hpp"

#include "BSG/Generator.h"

TEST_CASE("Basic construction") {
  std::string iniFile = "/home/leendert/git/BSG/tests/45Ca.ini";
  BSG::Generator gen;

  bool success = gen.InitializeTransitionFromFile(iniFile);

  REQUIRE(success == true);

  BSG::BetaParams bp = gen.GetBetaParams();

  REQUIRE(bp.Zi == 20);
  REQUIRE(bp.A == 45);
  REQUIRE(bp.Zf == 21);
}
