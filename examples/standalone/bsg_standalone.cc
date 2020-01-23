#include "BSG/Generator.h"

#include "CLI11.hpp"

int main(int argc, const char** argv) {
  std::string iniFilename;
  std::string configFilename;
  std::string outputName = "output";

  CLI::App app{"Beta Spectrum Generator standalone"};
  app.add_option("-i,--input", iniFilename, "INI input file for transition information")->required();
  app.add_option("-c,--config", configFilename, "INI config file for calculation information");
  app.add_option("-o,--output", outputName, "Name for file output. No extensions.");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    app.exit(e);
  }

  BSG::Generator gen = BSG::Generator(outputName);
  bool success = gen.Initialize(iniFilename, configFilename, argc, argv);
  if (success)
    gen.CalculateSpectrum();

  return 0;
}
