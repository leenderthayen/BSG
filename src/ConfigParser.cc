#include "NME/ConfigParser.h"

#include "CLI11.hpp"

namespace BSG {

  void parse(CLI::App& app, int argc, char** argv) {
    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      app.exit(e);
    }
  }

  void parseEmpty(CLI::App& app) {
    int argc = 1;
    char* argv[1] = {"coupling"};

    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      app.exit(e);
    }
  }

  CouplingConstants ParseCouplingConstants (std::string filename) {
    CouplingConstants couplingConstants;

    CLI::App app{"Coupling constants"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* coupling = app.add_subcommand("Coupling", "This is the coupling subcommand");
    coupling->add_option("--gV", couplingConstants.gV, "Vector coupling constant.");
    coupling->add_option("--gA", couplingConstants.gA, "Axial vector coupling constant.");
    coupling->add_option("--gM", couplingConstants.gM, "Weak magnetism coupling constant.");
    coupling->add_option("--gP", couplingConstants.gP, "Induced pseudoscalar coupling constant.");

    parseEmpty(app);

    return couplingConstants;
  }

  CorrectionOptions ParseCorrectionOptions (std::string filename) {
    CorrectionOptions correctionOptions;

    CLI::App app{"Correction options"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* correction = app.add_subcommand("Spectrum", "This is the spectrum subcommand");
    correction->add_option("-p,--phasespace", correctionOptions.phaseSpace, "Vector correction constant.");
    correction->add_option("-f,--fermi", correctionOptions.FermiFunction, "Axial vector correction constant.");
    correction->add_option("-l,--esfinitesize", correctionOptions.ESFiniteSize, "Weak magnetism correction constant.");
    correction->add_option("-C,--shapefactor", correctionOptions.shapeFactor, "");
    correction->add_option("-I,--isovector", correctionOptions.isovector, "");
    correction->add_option("-R,--relativistic", correctionOptions.relativistic, "");
    correction->add_option("-D,--esdeformation", correctionOptions.ESDeformation, "");
    correction->add_option("-U,--esfermi", correctionOptions.ESFermi, "");
    correction->add_option("-Q,--coulombrecoil", correctionOptions.CoulombRecoil, "");
    correction->add_option("-r,--radiative", correctionOptions.radiative, "");
    correction->add_option("-n,--kinematicrecoil", correctionOptions.kinRecoil, "");
    correction->add_option("-s,--screening", correctionOptions.atomicScreen, "");
    correction->add_option("-x,--exchange", correctionOptions.atomicExchange, "");
    correction->add_option("-m,--atomicmismatch", correctionOptions.atomicMismatch, "");

    parseEmpty(app);

    return correctionOptions;
  }

  SpectrumCalculationOptions ParseSpectrumCalculationOptions(filename) {
    SpectrumCalculationOptions spectrumCalcOptions;

    CLI::App app{"Spectrum calculation options"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* calc = app.add_subcommand("Calculation", "This is the calculation subcommand");
    calc->add_option("-b,--begin", spectrumCalcOptions.begin, "");
    calc->add_option("-e,--end", spectrumCalcOptions.end, "");
    calc->add_option("-s,--stepsize", spectrumCalcOptions.stepSize, "");
    calc->add_option("-n,--steps", spectrumCalcOptions.steps, "");
    calc->add_option("-v,--neutrino", spectrumCalcOptions.includeNeutrino, "");

    parseEmpty(app);

    return spectrumCalcOptions;
  }

  AdvancedOptions ParseAdvancedOptions(std::string filename) {
    AdvancedOptions advancedOptions;

    CLI::App app{"Advanced options"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* adv = app.add_subcommand("Advanced", "This is the advanced subcommand");
    adv->add_option("-c,--connect", spectrumCalcOptions.connectSPS, "");
    adv->add_option("-e,--esshape", spectrumCalcOptions.ESShape, "");
    adv->add_option("-n,--nsshape", spectrumCalcOptions.NSShape, "");
    adv->add_option("-v,--vold", spectrumCalcOptions.vold, "");
    adv->add_option("-V,--vnew", spectrumCalcOptions.vnew, "");

    parseEmpty(app);

    return advancedOptions;
  }

  AllowedMatrixElements ParseAllowedMatrixElements {
    AllowedMatrixElements allowedME;

    CLI::App app{"Allowed matrix elements"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* me = app.add_subcommand("MatrixElements", "This is the MatrixElements subcommand");
    me->add_option("-b,--weakmagnetism", allowedME.bAc, "");
    me->add_option("-d,--inducedtensor", allowedME.dAc, "");
    me->add_option("-l,--lambda", allowedME.lambda, "");

    parseEmpty(app);

    return allowedME;
  }

  ConfigOptions ParseConfigFile(std::string filename) {
    ConfigOptions configOptions;

    configOptions.couplingConstants = ParseCouplingConstants(filename);
    configOptions.correctionOptions = ParseCorrectionOptions(filename);
    configOptions.spectrumCalcOptions = ParseSpectrumCalculationOptions(filename);
    configOptions.advancedOptions = ParseAdvancedOptions(filename);
    configOptions.allowedME = ParseAllowedMatrixElements(filename);

    return configOptions;
  }
}//end of NME namespace
