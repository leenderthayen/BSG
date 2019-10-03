#include "catch.hpp"

#include "BSG/ConfigParser.h"
#include "PDS/Core/ParticleDefinition.h"
#include "PDS/Core/Nucleus.h"

TEST_CASE("Test general help") {
  const int argc = 2;
  const char* argv[argc] = {"test", "-h"};
  BSG::ParseConfigOptions("", argc, argv);
}

TEST_CASE("Test coupling help") {
  const int argc = 3;
  SECTION("lowercase") {
    const char* argv[argc] = {"test", "coupling", "-h"};
    BSG::ParseConfigOptions("", argc, argv);
  }
  SECTION("UPPERCASE") {
    const char* argv[argc] = {"test", "COUPLING", "-h"};
    BSG::ParseConfigOptions("", argc, argv);
  }
}
