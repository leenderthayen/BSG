#include "NME/ConfigParser.h"

#include "CLI11.hpp"

namespace NME {

  void parse(CLI::App& app) {
    int argc = 1;
    char* argv[1] = {"coupling"};

    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      app.exit(e);
    }
  }

  void ParseROBTDMethod(std::string filename, ConfigOptions& configOptions) {
    CLI::App app{"Coupling constants"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* comp = app.add_subcommand("Computational", "This is the computational subcommand");
    comp->add_option("--ROBTDMethod", configOptions.robtdMethod, "");

    parse(app);
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

    parse(app);

    return couplingConstants;
  }

  WSPotentialOptions ParseWSPotentialOptions (std::string filename) {
    WSPotentialOptions wsPotentialOptions;

    CLI::App app{"Woods-Saxon Potential Options"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* comp = app.add_subcommand("Computational", "This is the computational subcommand");
    CLI::App* ws = comp->add_subcommand("Woods-Saxon", "Woods-Saxon potential options");
    ws->add_option("--SurfaceThickness", wsPotentialOptions.surfaceThickness, "Surface thickness in fm");
    ws->add_option("--Vneutron", wsPotentialOptions.vNeutron, "");
    ws->add_option("--Vproton", wsPotentialOptions.vProton, "");
    ws->add_option("--Xneutron", wsPotentialOptions.xNeutron, "");
    ws->add_option("--Xproton", wsPotentialOptions.xProton, "");
    ws->add_option("--V0Sneutron", wsPotentialOptions.v0sNeutron, "");
    ws->add_option("--V0Sproton", wsPotentialOptions.v0sProton, "");

    parse(app);

    return wsPotentialOptions;
  }

  ESPOptions ParseESPOptions (std::string filename) {
    ESPOptions espOptions;
    CLI::App app{"Extremely Single Particle Options"};
    app.allow_config_extras(true);
    app.set_config("-c", filename);
    CLI::App* comp = app.add_subcommand("Computational", "This is the computational subcommand");
    CLI::App* esp = comp->add_subcommand("ESP", "Extremely Single Particle options");
    esp->add_option("--Potential", espOptions.potential, "");
    esp->add_option("--EnergyMargin", espOptions.energyMargin, "");
    esp->add_option("--ForceSpin", espOptions.forceSpin, "");
    esp->add_option("--ReversedGhallagher", espOptions.reversedGhallagher, "");
    esp->add_option("--OverrideSPCoupling", espOptions.overrideSPCoupling, "");

    parse(app);

    return espOptions;
  }

  ConfigOptions ParseConfigFile(std::string filename) {
    ConfigOptions configOptions;

    ParseROBTDMethod(filename, configOptions);

    configOptions.couplingConstants = ParseCouplingConstants(filename);
    configOptions.wsPotentialOptions = ParseWSPotentialOptions(filename);
    configOptions.espOptions = ParseESPOptions(filename);

    return configOptions;
  }
}//end of NME namespace
