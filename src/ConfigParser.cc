#include "BSG/ConfigParser.h"
#include "NHL/ParticleParser.h"

namespace BSG {

  void parse(CLI::App& app, int argc = 0, const char** argv = nullptr) {
    if (argc == 0) {
      argc = 1;
      const char* _argv[1] = {"coupling"};
      argv = _argv;
    }
    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      app.exit(e);
    }
  }

  void SetTransitionOptions(CLI::App& app, TransitionOptions& transitionOptions) {
    CLI::App* trans = app.add_subcommand("Transition", "This is the transition subcommand")->ignore_case();
    trans->add_option("-b,--betatype", transitionOptions.betaType, "")->ignore_case();
    trans->add_option("-d,--decaytype", transitionOptions.decayType, "")->ignore_case();
    trans->add_option("-m,--mixing", transitionOptions.mixingRatio, "")->ignore_case();
    trans->add_option("-Q,--qvalue", transitionOptions.QValue, "")->ignore_case();
    trans->add_option("-a,--atomicenergydeficit", transitionOptions.atomicEnergyDeficit, "")->ignore_case();
    trans->add_option("-p,--partialhalflife", transitionOptions.partialHalflife, "")->ignore_case();
    trans->add_option("-l,--logft", transitionOptions.logft, "")->ignore_case();
  }

  void SetCouplingOptions(CLI::App& app, CouplingConstants& couplingConstants) {
    CLI::App* coupling = app.add_subcommand("Coupling", "This is the coupling subcommand")->ignore_case();
    coupling->add_option("--gV", couplingConstants.gV, "Vector coupling constant.")->ignore_case();
    coupling->add_option("--gA", couplingConstants.gA, "Axial vector coupling constant.")->ignore_case();
    coupling->add_option("--gM", couplingConstants.gM, "Weak magnetism coupling constant.")->ignore_case();
    coupling->add_option("--gP", couplingConstants.gP, "Induced pseudoscalar coupling constant.")->ignore_case();
  }

  void SetCorrectionOptions(CLI::App& app, CorrectionOptions& correctionOptions) {
    CLI::App* correction = app.add_subcommand("Spectrum", "This is the spectrum subcommand")->ignore_case();
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
  }

  void SetSpectrumCalculationOptions(CLI::App& app, SpectrumCalculationOptions& spectrumCalcOptions) {
    CLI::App* calc = app.add_subcommand("Calculation", "This is the calculation subcommand")->ignore_case();
    calc->add_option("-b,--begin", spectrumCalcOptions.begin, "")->ignore_case();
    calc->add_option("-e,--end", spectrumCalcOptions.end, "")->ignore_case();
    calc->add_option("-s,--stepsize", spectrumCalcOptions.stepSize, "")->ignore_case();
    calc->add_option("-n,--steps", spectrumCalcOptions.steps, "")->ignore_case();
    calc->add_option("-v,--neutrino", spectrumCalcOptions.includeNeutrino, "")->ignore_case();
  }

  void SetAdvancedOptions(CLI::App& app, AdvancedOptions& advancedOptions) {
    CLI::App* adv = app.add_subcommand("Advanced", "This is the advanced subcommand")->ignore_case();
    adv->add_option("-c,--connect", advancedOptions.connectSPS, "")->ignore_case();
    adv->add_option("-e,--esshape", advancedOptions.ESShape, "")->ignore_case();
    adv->add_option("-n,--nsshape", advancedOptions.NSShape, "")->ignore_case();
    adv->add_option("-v,--vold", advancedOptions.vold, "");
    adv->add_option("-V,--vnew", advancedOptions.vnew, "");
  }

  void SetAllowedMatrixElementOptions(CLI::App& app, AllowedMatrixElements& allowedME) {
    CLI::App* me = app.add_subcommand("MatrixElements", "This is the MatrixElements subcommand")->ignore_case();
    me->add_option("-b,--weakmagnetism", allowedME.bAc, "")->ignore_case();
    me->add_option("-d,--inducedtensor", allowedME.dAc, "")->ignore_case();
    me->add_option("-l,--lambda", allowedME.lambda, "")->ignore_case();
  }

  TransitionOptions ParseTransitionOptions(std::string filename) {
    TransitionOptions transitionOptions;

    transitionOptions.initNucleus = NHL::ParseNucleus(filename, "Mother");
    transitionOptions.finalNucleus = NHL::ParseNucleus(filename, "Daughter");

    CLI::App app{"Beta Spectrum Generator Transition Options"};
    app.allow_extras(true);
    app.allow_config_extras(true);
    app.set_config("-c", filename);

    SetTransitionOptions(app, transitionOptions);

    parse(app);

    return transitionOptions;
  }

  ConfigOptions ParseConfigOptions(std::string filename, int argc, const char** argv) {
    ConfigOptions configOptions;

    CLI::App app{"Beta Spectrum Generator Configuration Options"};
    app.allow_extras(true);
    app.allow_config_extras(true);
    app.set_config("-c", filename);

    SetCouplingOptions(app, configOptions.couplingConstants);
    SetCorrectionOptions(app, configOptions.correctionOptions);
    SetSpectrumCalculationOptions(app, configOptions.spectrumCalcOptions);
    SetAdvancedOptions(app, configOptions.advancedOptions);
    SetAllowedMatrixElementOptions(app, configOptions.allowedME);

    parse(app, argc, argv);

    return configOptions;
  }

  bool TransitionSanityCheck(TransitionOptions& tr) {
    return true;
  }
}//end of BSG namespace
