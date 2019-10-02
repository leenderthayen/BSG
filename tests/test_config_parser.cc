#include "catch.hpp"

#include "BSG/ConfigParser.h"

TEST_CASE("Test general help") {
  const int argc = 2;
  const char* argv[argc] = {"test", "-h"};
  BSG::ParseOptions("", argc, argv);
}

TEST_CASE("Test coupling help") {
  const int argc = 3;
  SECTION("lowercase") {
    const char* argv[argc] = {"test", "coupling", "-h"};
    BSG::ParseOptions("", argc, argv);
  }
  SECTION("UPPERCASE") {
    const char* argv[argc] = {"test", "COUPLING", "-h"};
    BSG::ParseOptions("", argc, argv);
  }
}
